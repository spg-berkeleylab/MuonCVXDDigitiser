#include "ChargeClustersBuilder.h"

/* ****************************************************************************

    Find-Union Algorithm

   ************************************************************************* */

GridPartitionedSet::GridPartitionedSet(int x_s, int y_s) :
    x_size(x_s),
    y_size(y_s),
    data(x_size * y_size, 0)
{
    reset();
}

void GridPartitionedSet::reset()
{
    for (size_t k = 0; k < data.size(); k++) data[k] = k;
}

/*
 * Remove internal references in the grid, each cell contains the cluster ID
 */
void GridPartitionedSet::collapse()
{
    for (int h = 0; h < x_size; h++)
    {
        for (int k = 0; k < y_size; k++)
        {
            if (data[k] >= 0) find(h, k);
        }
    }
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

            _gridSet.reset();

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
            
            _gridSet.collapse();
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



