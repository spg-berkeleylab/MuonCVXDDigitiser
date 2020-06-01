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
                         int layer,
                         int ladder,
                         float tclick,
                         float wsize,
                         float thickness,
                         double tanLorentzAngleX,
                         double tanLorentzAngleY,
                         double cutOnDeltaRays,
                         double diffusionCoefficient,
                         double electronsPerKeV,
                         double segmentLength,
                         double energyLoss,
                         SurfaceMap* s_map);
    virtual ~DetElemSlidingWindow();
    void move_forward();
    PixelDigiMatrix getPixels();

private:
    void StoreSignalPoints(SimTrackerHit* hit);

    float curr_time;
    float time_click;
    float window_radius;

    HitTemporalIndexes& _htable;
    int _layer;
    int _ladder;
    float _layerThickness;
    float _layerHalfThickness;
    double _tanLorentzAngleX;
    double _tanLorentzAngleY;
    double _cutOnDeltaRays;
    double _diffusionCoefficient;
    double _electronsPerKeV;
    double _segmentLength;
    double _energyLoss;

    TimedSignalPointList signals;
    SurfaceMap* surf_map;
    G4UniversalFluctuation* _fluctuate;
};






#endif //DetElemSlidingWindow_h
