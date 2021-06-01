#ifndef HKBaseSensor_h
#define HKBaseSensor_h 1

#include "PixelDigiMatrix.h"
#include <tuple>
#include <unordered_map>

using std::vector;
using std::tuple;
using std::tie;
using std::unordered_map;

struct ClusterData
{
    int label;
    int pos;
};

static bool CmpClusterData(ClusterData c1, ClusterData c2) { return c1.label < c2.label; }

/* ****************************************************************************

    Find-Union Algorithm

   ************************************************************************* */

struct GridCoordinate
{
    int row;
    int col;
};

static inline bool operator==(GridCoordinate a, GridCoordinate b)
{
    return a.row == b.row && a.col == b.col;
}

using ClusterOfPixel = vector<GridCoordinate>;

tuple<int, int, int, int> GetBound(const ClusterOfPixel& cluster);

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
    ClusterOfPixel next();

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

/* ****************************************************************************

    Cluster Heap

   ************************************************************************* */

struct CounterItem
{
    int   left;
    int   head;
//    float time;
};

using CounterTable = unordered_map<int, CounterItem>;

struct ChargeItem
{
    float charge;
    int   next;
    int   cid;
};

using ChargeTable = unordered_map<int, ChargeItem>;

struct ChargePoint
{
    int row;
    int col;
    float charge;
};

struct BufferedCluster
{
    vector<ChargePoint> pixels;
    float time;
};

class ClusterHeap
{
public:
    ClusterHeap(int b_size);
    virtual ~ClusterHeap();
    void AddCluster(ClusterOfPixel& cluster);
    void UpdatePixel(int pos_x, int pos_y, PixelData pix);
    vector<BufferedCluster> PopClusters();
private:
    int hash_cnt;
    int bunch_size;
    ChargeTable  charge_table;
    CounterTable counter_table;
    vector<int>  ready_to_pop;
};

/* ****************************************************************************

    Hoshen-Kopelman sensor

   ************************************************************************* */

class HKBaseSensor : public PixelDigiMatrix
{
public:
    HKBaseSensor(int layer,
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
                          float t_step);
    virtual ~HKBaseSensor() {}

    void buildHits(SegmentDigiHitList& output) override;

protected:
    virtual bool pixelOn(int seg_x, int seg_y, int pos_x, int pos_y);
    virtual ClusterOfPixel processCluster(const ClusterOfPixel& in) { return in; };

    GridPartitionedSet  _gridSet;
    vector<ClusterHeap> heap_table;
};

#endif //HKBaseSensor_h

