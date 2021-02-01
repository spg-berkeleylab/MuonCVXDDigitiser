#include "HKOSensor.h"
#include <UTIL/BitField64.h>
#include <UTIL/LCTrackerConf.h>
#include <UTIL/ILDConf.h>

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

vector<GridCoordinate> GridPartitionedSet::next()
{
    c_curr = c_next;
    if (c_curr == c_buffer.size())
    {
        vector<GridCoordinate> empty {};
        return empty;
    }

    while (c_next < c_buffer.size() && c_buffer[c_next].label == c_buffer[c_curr].label) c_next++;

    int res_size = c_next - c_curr;
    vector<GridCoordinate> result(res_size, { 0, 0 });

    for (int k = 0; k < res_size; k++)
    {
        result[k] = coordinate(c_buffer[c_curr + k].pos);
    }
    return result;
}

/* ****************************************************************************

    Hoshen-Kopelman-Otsu sensor

   ************************************************************************* */

HKOSensor::HKOSensor(int layer,
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
                    barrel_id),
    _q_level(q_level),
    _gridSet(s_rows, s_colums)
{}

void HKOSensor::buildHits(SegmentDigiHitList& output)
{
    BitField64 bf_encoder(this->GetCellIDFormatStr());
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
                        _gridSet.merge(i - 1, j, i, j);
                    }
                    else if (up_is_above && !left_is_above)
                    {
                        _gridSet.merge(i, j - 1, i, j);
                    }
                }
            }
            
            _gridSet.close();

            vector<GridCoordinate> c_item = _gridSet.next();
            while (c_item.size() > 0)
            {
                // Very simple implemetation: center-of-gravity
                float tot_charge = 0;
                float x_acc = 0;
                float y_acc = 0;
                for (GridCoordinate p_coord : c_item)
                {
                    int global_x = SensorRowToLadderRow(h, p_coord.x);
                    int global_y = SensorColToLadderCol(k, p_coord.y);
                    float tmpc = GetPixel(global_x, global_y).charge;
                    double pos_x = PixelRowToX(global_x);
                    double pos_y = PixelColToY(global_y);

                    tot_charge += tmpc;
                    x_acc += tmpc * pos_x;
                    y_acc += tmpc * pos_y;
                }

                //Sensor segments ordered row first
                bf_encoder[LCTrackerCellID::sensor()] = h * this->GetSegNumX() + k;

                SegmentDigiHit digiHit = {
                    x_acc / tot_charge,
                    y_acc / tot_charge,
                    tot_charge,
                    0,                    // TODO missing sampling time
                    bf_encoder.lowWord(),
                    h, k
                };
                output.push_back(digiHit);

                c_item = _gridSet.next();
            }
        }
    }
}

float HKOSensor::getThreshold(int segid_x, int segid_y)
{
    // Otsu algorithm
    //   https://en.wikipedia.org/wiki/Otsu%27s_method
    //   http://www.labbookpages.co.uk/software/imgProc/otsuThreshold.html

    // Quantization based on max charge value over the ladder
    float max_chrg = GetMaxCharge();
    if (max_chrg == 0) return 0;

    float chrg_step = max_chrg / _q_level;

    // histogram
    vector<int> histo { _q_level, 0 };

    for (int h = 0; h < GetSensorRows(); h++)
    {
        for (int k = 0; k < GetSensorCols(); k++)
        {
            float tmpchrg = GetPixel(segid_x, segid_y, h, k).charge;

            int slot = int(floorf(tmpchrg / chrg_step));
            histo[slot]++;
        }
    }

    // Otsu algorithm
    int total_pixels = GetSensorRows() * GetSensorCols();
    float w_sum = 0;
    for (int j = 0; j < _q_level; j++) w_sum += j * histo[j];

    int threshold = 0;
    float varMax = 0;
    float sumB = 0;
    int wB = 0;
    int wF = 0;

    for (int j = 0; j < _q_level; j++)
    {
        wB += histo[j];
        if (wB == 0) continue;

        wF = total_pixels - wB;
        if (wF == 0) break;

        sumB += j * histo[j];

        float mB = sumB / wB;
        float mF = (w_sum - sumB) / wF;

        // Calculate Between Class Variance
        float varBetween = (float) wB * (float) wF * (mB - mF) * (mB - mF);

        // Check if new maximum found
        if (varBetween > varMax) {
            varMax = varBetween;
            threshold = j;
        }
    }

    return chrg_step * (threshold + 1);
}

bool HKOSensor::aboveThreshold(float charge, int seg_x, int seg_y, int pos_x, int pos_y)
{
    if (pos_x < 0 || pos_x >= this->GetSensorRows()) return false;
    if (pos_y < 0 || pos_y >= this->GetSensorCols()) return false;

    return GetPixel(seg_x, seg_y, pos_x, pos_y).charge > charge;
}




