#ifndef HKBaseSensor_h
#define HKBaseSensor_h 1

#include "PixelDigiMatrix.h"
#include <tuple>
#include <unordered_map>

using std::vector;
using std::tuple;
using std::tie;
using std::unordered_map;

/* ****************************************************************************

    Find-Union Algorithm

   ************************************************************************* */

using ClusterOfPixel = vector<LinearPosition>;

tuple<int, int, int, int> GetBound(const ClusterOfPixel& cluster, GridPosition locate);

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
    vector<ClusterOfPixel> get_clusters();

private:
    int rows;
    int columns;
    unordered_map<int, ClusterOfPixel> c_table;

    GridPosition locate;
    vector<int> data;
};

/* ****************************************************************************

    Cluster Heap

   ************************************************************************* */
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

struct ClusterItem
{
    BufferedCluster buffer;
    int size;
};

using ClusterTable = unordered_map<int, ClusterItem>;

using ReferenceTable = unordered_map<LinearPosition, int>;

class ClusterHeap
{
public:
    ClusterHeap(int rows, int cols);
    virtual ~ClusterHeap();
    void AddCluster(ClusterOfPixel& cluster);
    void SetupPixel(int pos_x, int pos_y, PixelData pix);
    vector<BufferedCluster> PopClusters();
    void SetLabel(string dlabel) { debug_label = dlabel; }

private:
    int hash_cnt;
    GridPosition locate;
    string debug_label;
    ClusterTable cluster_table;
    ReferenceTable ref_table;
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
    GridPartitionedSet  _gridSet;
    vector<ClusterHeap> heap_table;
};

#endif //HKBaseSensor_h

