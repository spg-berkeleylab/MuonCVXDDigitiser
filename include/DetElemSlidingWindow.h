#ifndef DetElemSlidingWindow_h
#define DetElemSlidingWindow_h 1

#include <list>

#include "HitTemporalIndexes.h"
#include "PixelDigiMatrix.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "G4UniversalFluctuation.h"

using dd4hep::rec::SurfaceMap;

struct TimedSignalPoint
{
    double x;
    double y;
    double sigmaX;
    double sigmaY;
    double charge;
    float time;
};

typedef std::list<TimedSignalPoint> TimedSignalPointList;

class DetElemSlidingWindow
{
public:
    DetElemSlidingWindow(HitTemporalIndexes& htable,
                         PixelDigiMatrix& sensor,
                         float tclick,
                         float wsize,
                         double tanLorentzAngleX,
                         double tanLorentzAngleY,
                         double cutOnDeltaRays,
                         double diffusionCoefficient,
                         double electronsPerKeV,
                         double segmentLength,
                         double energyLoss,
                         double widthOfCluster,
                         SurfaceMap* s_map);
    virtual ~DetElemSlidingWindow();
    bool move_forward();

private:
    void StoreSignalPoints(SimTrackerHit* hit);
    void UpdatePixels();

    float curr_time;
    float time_click;
    float window_radius;

    HitTemporalIndexes& _htable;
    PixelDigiMatrix& _sensor;
    double _tanLorentzAngleX;
    double _tanLorentzAngleY;
    double _cutOnDeltaRays;
    double _diffusionCoefficient;
    double _electronsPerKeV;
    double _segmentLength;
    double _energyLoss;
    double _widthOfCluster;
    TimedSignalPointList signals;
    SurfaceMap* surf_map;
    G4UniversalFluctuation* _fluctuate;
};






#endif //DetElemSlidingWindow_h
