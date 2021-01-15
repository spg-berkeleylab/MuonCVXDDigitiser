#ifndef ChargeClustersBuilder_h
#define ChargeClustersBuilder_h 1

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
    GridPartitionedSet(int x_s, int y_s);
    virtual ~GridPartitionedSet() {}
    int find(int x, int y);
    void merge(int x1, int y1, int x2, int y2);
    void init();
    void close();
    void invalidate(int x, int y);
    vector<GridCoordinate> next();

private:
    inline int index(int x, int y) { return x * x_size + y; }
    inline GridCoordinate coordinate(int p) { return { p / x_size, p % x_size }; }

    int x_size;
    int y_size;
    int valid_cells;
    int c_curr;
    int c_next;
    vector<int> data;
    vector<ClusterData> c_buffer;
};

class ChargeClustersBuilder : public PixelDigiMatrix
{
public:
    ChargeClustersBuilder(int layer,
                          int ladder,
                          int xsegmentNumber,
                          int ysegmentNumber,
                          float ladderLength,
                          float ladderWidth,
                          float thickness,
                          double pixelSizeX,
                          double pixelSizeY,
                          string enc_str);
    virtual ~ChargeClustersBuilder() {}

    void buildHits(SegmentDigiHitList& output);

protected:
    virtual float getThreshold(int segid_x, int segid_y);
    virtual bool aboveThreshold(float charge, int seg_x, int seg_y, int pos_x, int pos_y);

    GridPartitionedSet _gridSet;

private:
    inline int sensor_posx(int seg_x, int pos_x) { return seg_x * this->GetSegSizeX() + pos_x; }
    inline int sensor_posy(int seg_y, int pos_y) { return seg_y * this->GetSegSizeY() + pos_y; }
};

#endif //ChargeClustersBuilder_h

