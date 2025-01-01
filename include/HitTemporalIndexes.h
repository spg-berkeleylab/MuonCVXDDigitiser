#ifndef HitTemporalIndexes_h
#define HitTemporalIndexes_h 1

#include <vector>
#include <unordered_map>
#include <queue>
#include <limits>

#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include "BitField64.h"

using std::priority_queue;
using std::unordered_map;
using edm4hep::SimTrackerHit;
using edm4hep::SimTrackerHitCollection;

class CmpTrackTime
{
public:
    CmpTrackTime(){}
    bool operator()(SimTrackerHit &alfa, SimTrackerHit &beta)
    {
        return alfa.getTime() > beta.getTime();

    }
};

typedef priority_queue<SimTrackerHit*, std::vector<SimTrackerHit*>, CmpTrackTime> hit_queue;


class HitTemporalIndexes
{
public:
    HitTemporalIndexes(const SimTrackerHitCollection STHcol);
    virtual ~HitTemporalIndexes();
    SimTrackerHit* CurrentHit(int layer, int ladder);
    void DisposeHit(int layer, int ladder);
    int GetHitNumber(int layer, int ladder);
    float GetMinTime();
    float GetMinTime(int layer, int ladder);

    static float MAXTIME;

private:
    inline int GetKey(int layer, int ladder);

    BitField64 cellid_decoder;
    unordered_map<int, hit_queue*> htable;
};

#endif //HitTemporalIndexes_h

