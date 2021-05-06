#ifndef ShapeProcessingSensor_h
#define ShapeProcessingSensor_h 1

#include "HKBaseSensor.h"
#include <vector>

using std::vector;

class ShapeProcessingSensor : public HKBaseSensor
{
public:
    ShapeProcessingSensor(int layer,
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
                          int q_level);
    virtual ~ShapeProcessingSensor() {}

protected:
    ClusterOfPixel processCluster(const ClusterOfPixel& in) override;
    vector<GridCoordinate> GetContour(const ClusterOfPixel& spot);

private:
    GridCoordinate GetNextPoint(GridCoordinate c, GridCoordinate p);
};

#endif //ShapeProcessingSensor_h