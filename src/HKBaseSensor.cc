#include "HKBaseSensor.h"
#include <UTIL/BitField64.h>
#include <UTIL/LCTrackerConf.h>
#include <UTIL/ILDConf.h>

#include "streamlog/streamlog.h"

#include <math.h>
#include <algorithm>

using std::sort;
using UTIL::BitField64;
using lcio::LCTrackerCellID;
using lcio::ILDDetID;

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
    c_curr = c_next;
    if (c_curr == c_buffer.size())
    {
        ClusterOfPixel empty { {}, 0, 0, 0, 0 };
        return empty;
    }

    while (c_next < c_buffer.size() && c_buffer[c_next].label == c_buffer[c_curr].label) c_next++;

    int res_size = c_next - c_curr;
    ClusterOfPixel result { { res_size, { 0, 0 } }, rows + 1, -1, columns + 1, -1 };

    for (int k = 0; k < res_size; k++)
    {
        GridCoordinate coord = coordinate(c_buffer[c_curr + k].pos);
        result.pix[k] = coord;
        if (coord.row < result.row_min) result.row_min = coord.row;
        if (coord.row > result.row_max) result.row_max = coord.row;
        if (coord.col < result.col_min) result.col_min = coord.col;
        if (coord.col > result.col_max) result.col_max = coord.col;
    }
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
                     float s_level,
                     int q_level) :
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
                    s_level,
                    q_level),
    _gridSet(s_rows, s_colums)
{}

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
            float charge_thr = getThreshold(h, k);
            if (charge_thr == 0) continue;

            //Hoshen-Kopelman Algorithm:
            //  https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html

            _gridSet.init();

            for (int i = 0; i < this->GetSensorRows(); i++)
            {
                for (int j = 0; j < this->GetSensorCols(); j++)
                {
                    if (!aboveThreshold(charge_thr, h, k, i, j))
                    {
                        _gridSet.invalidate(i, j);
                        continue;
                    }

                    bool up_is_above = aboveThreshold(charge_thr, h, k, i - 1, j);
                    bool left_is_above = aboveThreshold(charge_thr, h, k, i, j - 1);

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

            for (ClusterOfPixel c_item = processCluster(_gridSet.next());
                 c_item.pix.size() > 0;
                 c_item = processCluster(_gridSet.next()))
            {
                // Very simple implementation: geometric mean
                float n_acc = 0;
                float x_acc = 0;
                float y_acc = 0;
                float t_acc = 0;
                float tot_charge = 0;
                for (GridCoordinate p_coord : c_item.pix)
                {
                    int global_x = SensorRowToLadderRow(h, p_coord.row);
                    int global_y = SensorColToLadderCol(k, p_coord.col);

                    n_acc += 1;
                    x_acc += PixelRowToX(global_x);
                    y_acc += PixelColToY(global_y);
                    PixelData p_data = GetPixel(global_x, global_y);
                    t_acc += p_data.time;
                    tot_charge += p_data.charge;
                }

                //Sensor segments ordered row first
                bf_encoder[LCTrackerCellID::sensor()] = h * this->GetSegNumY() + k;

                SegmentDigiHit digiHit = {
                    x_acc / n_acc,
                    y_acc / n_acc,
                    tot_charge,
                    t_acc / n_acc,
                    bf_encoder.lowWord(),
                    h, k
                };
                output.push_back(digiHit);
            }
        }
    }
}

float HKBaseSensor::getThreshold(int segid_x, int segid_y)
{
	// TODO implement noisy pixels and threshold dispersion effects
    return _thr_level;
}

bool HKBaseSensor::aboveThreshold(float charge, int seg_x, int seg_y, int pos_x, int pos_y)
{
    if (pos_x < 0 || pos_x >= this->GetSensorRows()) return false;
    if (pos_y < 0 || pos_y >= this->GetSensorCols()) return false;

    return GetPixel(seg_x, seg_y, pos_x, pos_y).charge > charge;
}




