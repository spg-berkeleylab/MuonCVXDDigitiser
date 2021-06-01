#include "HKBaseSensor.h"
#include <UTIL/BitField64.h>
#include <UTIL/LCTrackerConf.h>
#include <UTIL/ILDConf.h>

#include "streamlog/streamlog.h"

#include <math.h>
#include <algorithm>
#include <limits>

using std::sort;
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
    valid_cells(0),
    c_curr(0),
    c_next(0),
    data(rows * columns, 0),
    c_buffer(0)
{}

void GridPartitionedSet::init()
{
    for (size_t k = 0; k < data.size(); k++) data[k] = k;
    valid_cells = data.size();
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
            if (data[index(h, k)] >= 0) find(h, k);
        }
    }

    /*
     * Prepare the ordered set of clusters
     */
    c_buffer.resize(valid_cells);
    int c_id = 0;
    for (size_t k = 0; k < data.size(); k++)
    {
        if (data[k] >= 0)
        {
            c_buffer[c_id] = { data[k], k };
            c_id++;
        }
    }
    sort(c_buffer.begin(), c_buffer.end(), CmpClusterData);

    c_curr = -1;
    c_next = 0;
}

int GridPartitionedSet::find(int x, int y)
{
    int res = index(x, y);
    while (data[res] != res)
    {
        if (res < 0) return -1;
        res = data[res];
    }

    int pos = index(x, y);
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
    data[index(x, y)] = -1;
    valid_cells--;
}

ClusterOfPixel GridPartitionedSet::next()
{
    ClusterOfPixel result;

    c_curr = c_next;
    if (c_curr == c_buffer.size())
    {
        return result;
    }

    while (c_next < c_buffer.size() && c_buffer[c_next].label == c_buffer[c_curr].label) c_next++;

    int res_size = c_next - c_curr;
    result.assign(res_size, { 0, 0 });

    for (int k = 0; k < res_size; k++)
    {
        GridCoordinate coord = coordinate(c_buffer[c_curr + k].pos);
        result[k] = coord;
    }
    return result;
}

tuple<int, int, int, int> GetBound(const ClusterOfPixel& cluster)
{
    int row_min = int_limits::max();
    int row_max = -1;
    int col_min = int_limits::max();
    int col_max = -1;
    for (auto item : cluster)
    {
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

ClusterHeap::ClusterHeap(int b_size) :
    hash_cnt(0),
    bunch_size(b_size),
    charge_table(),
    counter_table(),
    ready_to_pop()
{}

ClusterHeap::~ClusterHeap()
{}

void ClusterHeap::AddCluster(ClusterOfPixel& cluster)
{
    int prev_pos = -1;
    int curr_pos = -1;
    for (GridCoordinate p_coord : cluster)
    {
        curr_pos = p_coord.row * bunch_size + p_coord.col;
        ChargeItem item { 0, prev_pos, hash_cnt };
        charge_table.emplace(curr_pos, item);   // TODO check for pix duplicated
        prev_pos = curr_pos;
    }

    CounterItem cnt_item { cluster.size(), curr_pos };
    counter_table.emplace(hash_cnt, cnt_item);

    hash_cnt++;  // TODO possible overflow
}

void ClusterHeap::UpdatePixel(int pos_x, int pos_y, PixelData pix)
{
    if (pix.status == PixelStatus::ready)
    {
        int pos = pos_x * bunch_size + pos_y;
        auto c_item = charge_table.find(pos);
        if (c_item != charge_table.end() and (c_item->second).charge == 0)
        {
            charge_table[pos].charge = pix.charge;
            int cluster_id = (c_item->second).cid;
            auto cnt_item = counter_table.find(cluster_id);
            if (cnt_item != counter_table.end())
            {
                counter_table[cluster_id].left -= 1;
                if (counter_table[cluster_id].left == 0)
                {
                    ready_to_pop.push_back(cluster_id);
                }
            }
        }
    }
}

vector<BufferedCluster> ClusterHeap::PopClusters()
{
    vector<BufferedCluster> result;
    for (auto cluster_id : ready_to_pop)
    {
        BufferedCluster c_points;

        auto cnt_item =  counter_table.find(cluster_id);
        if (cnt_item != counter_table.end())
        {
            int curr_pos = (cnt_item->second).head;
            do
            {
                auto ch_item = charge_table.find(curr_pos);
                if (ch_item != charge_table.end())
                {
                    ChargePoint c_pix
                    {
                        curr_pos / bunch_size,
                        curr_pos % bunch_size,
                        (ch_item->second).charge
                    };
                    
                    c_points.pixels.push_back(c_pix);
                }
                else
                {
                    curr_pos = -1;
                }
            }
            while(curr_pos != -1);
        }
        result.push_back(c_points);
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
    heap_table(0, { 0 })
{
    if (GetStatus() == MatrixStatus::ok)
    {
        heap_table.resize(GetSegNumX() * GetSegNumY(), { GetSensorCols() });
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
            /* ****************************************************************
               Hoshen-Kopelman Algorithm:
               https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html
               ************************************************************** */

            _gridSet.init();

            for (int i = 0; i < this->GetSensorRows(); i++)
            {
                for (int j = 0; j < this->GetSensorCols(); j++)
                {
                    if (!pixelOn(h, k, i, j))
                    {
                        _gridSet.invalidate(i, j);
                        continue;
                    }

                    bool up_is_above = pixelOn(h, k, i - 1, j);
                    bool left_is_above = pixelOn(h, k, i, j - 1);

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
            ClusterHeap& c_heap = heap_table[h * GetSegNumY() + k];
            for (ClusterOfPixel c_item = _gridSet.next();
                                c_item.size() > 0;
                                c_item = _gridSet.next())
            {
                c_heap.AddCluster(c_item);
            }

            for (int i = 0; i < this->GetSensorRows(); i++)
            {
                for (int j = 0; j < this->GetSensorCols(); j++)
                {
                    c_heap.UpdatePixel(i, j, GetPixel(h, k, i, j));
                }
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

                //Sensor segments ordered row first
                bf_encoder[LCTrackerCellID::sensor()] = h * this->GetSegNumY() + k;

                SegmentDigiHit digiHit = {
                    x_acc / c_item.pixels.size(),
                    y_acc / c_item.pixels.size(),
                    tot_charge,
                    c_item.time,
                    bf_encoder.lowWord(),
                    h, k
                };
                output.push_back(digiHit);
            }
        }
    }
}

bool HKBaseSensor::pixelOn(int seg_x, int seg_y, int pos_x, int pos_y)
{
    /*
    if (pos_x < 0 || pos_x >= this->GetSensorRows()) return false;
    if (pos_y < 0 || pos_y >= this->GetSensorCols()) return false;
    */
    return GetPixel(seg_x, seg_y, pos_x, pos_y).status == PixelStatus::start;
}




