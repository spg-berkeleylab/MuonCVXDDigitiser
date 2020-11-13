#include "ChargeClustersBuilder.h"

#include <algorithm>

using std::sort;

/* ****************************************************************************

    Find-Union Algorithm

   ************************************************************************* */

GridPartitionedSet::GridPartitionedSet(int x_s, int y_s) :
    x_size(x_s),
    y_size(y_s),
    valid_cells(0),
    c_curr(0),
    c_next(0),
    data(x_size * y_size, 0),
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
    for (int h = 0; h < x_size; h++)
    {
        for (int k = 0; k < y_size; k++)
        {
            if (data[k] >= 0) find(h, k);
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
            c_buffer[c_id] = { k, data[k] };
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

    if (pset1 < pset2)
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

    Cluster builder

   ************************************************************************* */

ChargeClustersBuilder::ChargeClustersBuilder(PixelDigiMatrix& sensor) :
    _sensor(sensor),
    _gridSet(_sensor.GetSegSizeX(), _sensor.GetSegSizeY())
{}

void ChargeClustersBuilder::buildHits(TrackerHitList& output)
{
    for (int h = 0; h < _sensor.GetSegNumX(); h++)
    {
        for (int k = 0; k < _sensor.GetSegNumY(); k++)
        {
            float charge_thr = getThreshold(h, k);

            //Hoshen-Kopelman Algorithm:
            //  https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html

            _gridSet.init();

            for (int i = 0; i < _sensor.GetSegSizeX(); i++)
            {
                for (int j = 0; j < _sensor.GetSegSizeY(); j++)
                {
                    if (!aboveThreshold(charge_thr, h, k, i, j))
                    {
                        _gridSet.invalidate(i, j);
                        continue;
                    }

                    bool up_is_above = i == 0 ? false : aboveThreshold(charge_thr, h, k, i - 1, j);
                    bool left_is_above = j == 0 ? false : aboveThreshold(charge_thr, h, k, i, j - 1);

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

            // TODO calculate tracker hits
        }
    }
}

float ChargeClustersBuilder::getThreshold(int segid_x, int segid_y)
{
    // TODO implement Otsu algorithm
    //   https://en.wikipedia.org/wiki/Otsu%27s_method
    return 0;
}

bool ChargeClustersBuilder::aboveThreshold(float charge, int seg_x, int seg_y, int pos_x, int pos_y)
{
    int global_x = seg_x * _sensor.GetSegSizeX() + pos_x;
    int global_y = seg_y * _sensor.GetSegSizeY() + pos_y;
    return _sensor.GetPixel(global_x, global_y).charge > charge;
}



