#include "HKBaseSensor.h"
#include <UTIL/BitField64.h>
#include <UTIL/LCTrackerConf.h>
#include <UTIL/ILDConf.h>

#include "streamlog/streamlog.h"

#include <math.h>
#include <algorithm>
#include <limits>
#include <sstream>

using std::sort;
using std::stringstream;
using UTIL::BitField64;
using lcio::LCTrackerCellID;
using lcio::ILDDetID;
using int_limits = std::numeric_limits<int>;

/* ****************************************************************************

    Find-Union Algorithm

   ************************************************************************* */

GridPartitionedSet::GridPartitionedSet(int n_row, int n_col) :
    rows(n_row),
    columns(n_col),
    c_table(),
    locate(n_row, n_col),
    data(rows * columns, 0)
{}

void GridPartitionedSet::init()
{
    for (size_t k = 0; k < data.size(); k++) data[k] = k;
    c_table.clear();
}

void GridPartitionedSet::close()
{
    /*
     * Remove internal references in the grid, each cell contains the cluster ID
     */
    for (int h = 0; h < rows; h++)
    {
        for (int k = 0; k < columns; k++)
        {
            if (data[locate(h, k)] >= 0) find(h, k);
        }
    }

    /*
     * Prepare the set of clusters
     */
    for (size_t k = 0; k < data.size(); k++)
    {
        auto c_item = c_table.find(data[k]);
        if (c_item != c_table.end())
        {
            (c_item->second).push_back(k);
        }
        else if (data[k] >= 0)
        {
            ClusterOfPixel p_list { k };
            c_table.emplace(data[k], p_list);
        }
    }
}

int GridPartitionedSet::find(int x, int y)
{
    int res = locate(x, y);
    while (data[res] != res)
    {
        if (res < 0) return -1;
        res = data[res];
    }

    int pos = locate(x, y);
    while (data[pos] != pos)
    {
        int curr = data[pos];
        data[pos] = res;
        pos = curr;
    }
    return res;
}

void GridPartitionedSet::merge(int x1, int y1, int x2, int y2)
{
    int pset1 = find(x1, y1);
    int pset2 = find(x2, y2);

    if (pset1 < 0 || pset2 < 0) return;

    if (pset1 >= pset2)
    {
        data[pset1] = pset2;
    }
    else
    {
        data[pset2] = pset1;
    }
}

void GridPartitionedSet::invalidate(int x, int y)
{
    data[locate(x, y)] = -1;
}

vector<ClusterOfPixel> GridPartitionedSet::get_clusters()
{
    vector<ClusterOfPixel> result;
    for (auto c_item : c_table)
    {
        result.push_back(std::move(c_item.second));
    }
    return result;
}

tuple<int, int, int, int> GetBound(const ClusterOfPixel& cluster, GridPosition locate)
{
    int row_min = int_limits::max();
    int row_max = -1;
    int col_min = int_limits::max();
    int col_max = -1;
    for (LinearPosition pos : cluster)
    {
        GridCoordinate item = locate(pos);
        if (item.row < row_min) row_min = item.row;
        if (item.row > row_max) row_max = item.row;
        if (item.col < col_min) col_min = item.col;
        if (item.col > col_max) col_max = item.col;
    }
    return std::make_tuple(row_min, row_max, col_min, col_max);
}

/* ****************************************************************************

    Cluster Heap

   ************************************************************************* */

ClusterHeap::ClusterHeap(int rows, int cols) :
    hash_cnt(0),
    locate(rows, cols),
    debug_label("Undefined"),
    cluster_table(),
    ref_table(),
    ready_to_pop()
{}

ClusterHeap::~ClusterHeap()
{}

void ClusterHeap::AddCluster(ClusterOfPixel& cluster)
{
    for (LinearPosition curr_pos : cluster)
    {
        auto ref_item = ref_table.find(curr_pos);
        if (ref_item == ref_table.end())
        {
            ref_table.emplace(curr_pos, hash_cnt);
        }
        else if (streamlog::out.write<streamlog::ERROR>())
#pragma omp critical
        {
            GridCoordinate gcoord = locate(curr_pos);
            streamlog::out() << "Cluster heap " << debug_label << ": conflict for pixel "
                << gcoord.row << ":" << gcoord.col << std::endl;
        }
    }

    ClusterItem cl_item { {}, cluster.size() };
    cluster_table.emplace(hash_cnt, cl_item);

    hash_cnt++;
    if (hash_cnt == int_limits::max())
    {
        hash_cnt = 0;
    }
}

void ClusterHeap::SetupPixel(int pos_x, int pos_y, PixelData pix)
{
    LinearPosition pos = locate(pos_x, pos_y);
    auto r_item = ref_table.find(pos);
    if (r_item != ref_table.end())
    {
        int cluster_id = r_item->second;
        ClusterItem& c_item = cluster_table[cluster_id];

        ChargePoint c_pix { pos_x, pos_y, pix.charge };
        c_item.buffer.pixels.push_back(c_pix);
        c_item.buffer.time = pix.time;

        if (c_item.buffer.pixels.size() == c_item.size)
        {
            ready_to_pop.push_back(cluster_id);
        }

        ref_table.erase(pos);
    }
    else if (streamlog::out.write<streamlog::ERROR>())
#pragma omp critical
    {
        streamlog::out() << "Cluster heap " << debug_label << ": undefined pixel "
            << pos_x << ":" << pos_y << std::endl;
    }
}

vector<BufferedCluster> ClusterHeap::PopClusters()
{
    vector<BufferedCluster> result;
    for (auto cluster_id : ready_to_pop)
    {
        result.push_back(cluster_table[cluster_id].buffer);
        cluster_table.erase(cluster_id);
    }

    ready_to_pop.clear();
    return result;
}

/* ****************************************************************************

    Hoshen-Kopelman sensor

   ************************************************************************* */

HKBaseSensor::HKBaseSensor(int layer,
                     int ladder,
                     int xsegmentNumber,
                     int ysegmentNumber,
                     float ladderLength,
                     float ladderWidth,
                     float thickness,
                     double pixelSizeX,
                     double pixelSizeY,
                     string enc_str,
                     int barrel_id,
                     double thr,
                     float fe_slope,
                     float starttime,
                     float t_step) :
    PixelDigiMatrix(layer,
                    ladder,
                    xsegmentNumber,
                    ysegmentNumber,
                    ladderLength,
                    ladderWidth,
                    thickness,
                    pixelSizeX,
                    pixelSizeY,
                    enc_str,
                    barrel_id,
                    thr,
                    fe_slope,
                    starttime,
                    t_step),
    _gridSet(s_rows, s_colums),
    heap_table(0, { 0, 0 })
{
    if (GetStatus() == MatrixStatus::ok)
    {
        heap_table.resize(GetSegNumX() * GetSegNumY(), { GetSensorRows(), GetSensorCols() });

        for (int h = 0; h < this->GetSegNumX(); h++)
        {
            for (int k = 0; k < this->GetSegNumY(); k++)
            {
                stringstream d_label;
                d_label << "[" << GetLayer() << ":" << GetLadder() << ":" << h << ":" << k << "]";
                heap_table[s_locate(h, k)].SetLabel(d_label.str());
            }
        }
    }
}

void HKBaseSensor::buildHits(SegmentDigiHitList& output)
{
    BitField64 bf_encoder { cellFmtStr };
    bf_encoder.reset();
    bf_encoder[LCTrackerCellID::subdet()] = _barrel_id;
    bf_encoder[LCTrackerCellID::side()] = ILDDetID::barrel;
    bf_encoder[LCTrackerCellID::layer()] = this->GetLayer();
    bf_encoder[LCTrackerCellID::module()] = this->GetLadder();

    for (int h = 0; h < this->GetSegNumX(); h++)
    {
        for (int k = 0; k < this->GetSegNumY(); k++)
        {

            //Sensor segments ordered row first
            LinearPosition sens_id = s_locate(h, k);
            bf_encoder[LCTrackerCellID::sensor()] = sens_id;

            ClusterHeap& c_heap = heap_table[sens_id];

            if (CheckStatusOnSensor(h, k, PixelStatus::start))
            {
                /* ****************************************************************
                   Hoshen-Kopelman Algorithm:
                   https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html
                   ************************************************************** */

                _gridSet.init();

                for (int i = 0; i < this->GetSensorRows(); i++)
                {
                    for (int j = 0; j < this->GetSensorCols(); j++)
                    {
                        if (!CheckStatus(h, k, i, j, PixelStatus::start))
                        {
                            _gridSet.invalidate(i, j);
                            continue;
                        }

                        bool up_is_above = CheckStatus(h, k, i - 1, j, PixelStatus::start);
                        bool left_is_above = CheckStatus(h, k, i, j - 1, PixelStatus::start);

                        if (up_is_above && left_is_above)
                        {
                            _gridSet.merge(i - 1, j, i, j - 1);
                            _gridSet.merge(i, j - 1, i, j);
                        }
                        else if (!up_is_above && left_is_above)
                        {
                            _gridSet.merge(i, j - 1, i, j);
                        }
                        else if (up_is_above && !left_is_above)
                        {
                            _gridSet.merge(i - 1, j, i, j);
                        }
                    }
                }

                _gridSet.close();

                /* ****************************************************************
                   Cluster buffering
                   ************************************************************** */

                for (ClusterOfPixel c_item : _gridSet.get_clusters())
                {
                    c_heap.AddCluster(c_item);
                }
            }

            for (auto p_item : GetPixelsFromSensor(h, k, PixelStatus::ready))
            {
                c_heap.SetupPixel(p_item.row, p_item.col, p_item.data);
            }

            for (BufferedCluster c_item : c_heap.PopClusters())
            {
                // Very simple implementation: geometric mean
                float x_acc = 0;
                float y_acc = 0;
                float tot_charge = 0;
                for (ChargePoint c_point : c_item.pixels)
                {
                    int global_x = SensorRowToLadderRow(h, c_point.row);
                    int global_y = SensorColToLadderCol(k, c_point.col);

                    x_acc += PixelRowToX(global_x);
                    y_acc += PixelColToY(global_y);

                    tot_charge += c_point.charge;
                }

                SegmentDigiHit digiHit = {
                    x_acc / c_item.pixels.size(),
                    y_acc / c_item.pixels.size(),
                    tot_charge,
                    c_item.time,
                    bf_encoder.lowWord()
                };
                output.push_back(digiHit);
            }
        }
    }
}
