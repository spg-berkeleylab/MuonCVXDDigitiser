#ifndef HKOSensor_h
#define HKOSensor_h 1

#include "PixelDigiMatrix.h"

using std::vector;

struct ClusterData
{
    int label;
    int pos;
};

static bool CmpClusterData(ClusterData c1, ClusterData c2) { return c1.label < c2.label; }

struct GridCoordinate
{
    int x;
    int y;
};

class GridPartitionedSet
{
public:
    GridPartitionedSet(int n_row, int n_col);
    virtual ~GridPartitionedSet() {}
    int find(int x, int y);
    void merge(int x1, int y1, int x2, int y2);
    void init();
    void close();
    void invalidate(int x, int y);
    vector<GridCoordinate> next();

private:
    inline int index(int row, int col) { return row * columns + col; }
    inline GridCoordinate coordinate(int p) { return { p / columns, p % columns }; }

    int rows;
    int columns;
    int valid_cells;
    int c_curr;
    int c_next;
    vector<int> data;
    vector<ClusterData> c_buffer;
};

class HKOSensor : public PixelDigiMatrix
{
public:
    HKOSensor(int layer,
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
                          int q_level);
    virtual ~HKOSensor() {}

    void buildHits(SegmentDigiHitList& output);

protected:
    virtual float getThreshold(int segid_x, int segid_y);
    virtual bool aboveThreshold(float charge, int seg_x, int seg_y, int pos_x, int pos_y);

    int _q_level;
    GridPartitionedSet _gridSet;
};

#endif //HKOSensor_h

