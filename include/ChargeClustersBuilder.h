#ifndef ChargeClustersBuilder_h
#define ChargeClustersBuilder_h 1

#include "PixelDigiMatrix.h"
#include <IMPL/TrackerHitPlaneImpl.h>

using std::vector;
using IMPL::TrackerHitPlaneImpl;

typedef vector<TrackerHitPlaneImpl> TrackerHitList;

struct ClusterData
{
    int label;
    int pos;
};

bool CmpClusterData(ClusterData c1, ClusterData c2) { return c1.label < c2.label; }

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

class ChargeClustersBuilder
{
public:
    ChargeClustersBuilder(PixelDigiMatrix& sensor);
    virtual ~ChargeClustersBuilder() {}

    virtual void buildHits(TrackerHitList& output);

protected:
    virtual float getThreshold(int segid_x, int segid_y);
    virtual bool aboveThreshold(float charge, int seg_x, int seg_y, int pos_x, int pos_y);

    PixelDigiMatrix& _sensor;
    GridPartitionedSet _gridSet;
};

#endif //ChargeClustersBuilder_h

