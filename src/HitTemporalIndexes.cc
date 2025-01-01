#include "HitTemporalIndexes.h"

using std::min;

HitTemporalIndexes::HitTemporalIndexes(const edm4hep::SimTrackerHitCollection STHcol):
    cellid_decoder("subdet:5,side:-2,layer:9,module:8,sensor:8"),
    htable()
{
    for (int i = 0; i < STHcol.size(); ++i)
    {
        SimTrackerHit simTrkHit = STHcol.at(i);
        cellid_decoder.setValue(simTrkHit.getCellID0());
        int layer = cellid_decoder["layer"];
        int ladder = cellid_decoder["module"];
        int tkey = GetKey(layer, ladder);

        auto item = htable.find(tkey);
        if (item == htable.end())
        {
            hit_queue* n_queue = new hit_queue();
            n_queue->push(&simTrkHit);
            htable.emplace(tkey, n_queue);
        }
        else
        {
            item->second->push(&simTrkHit);
        }
    }
}

HitTemporalIndexes::~HitTemporalIndexes()
{
    for (auto item : htable)
    {
        delete(item.second);
    }
}

SimTrackerHit HitTemporalIndexes::CurrentHit(int layer, int ladder)
{
    int tkey = GetKey(layer, ladder);
    auto item = htable.find(tkey);
    if (item != htable.end() && !item->second->empty())
    {
        return item->second->top();
    }
    return nullptr;
}

void HitTemporalIndexes::DisposeHit(int layer, int ladder)
{
    int tkey = GetKey(layer, ladder);
    auto item = htable.find(tkey);
    if (item != htable.end() && !item->second->empty())
    {
        item->second->pop();
    }
}

int HitTemporalIndexes::GetHitNumber(int layer, int ladder)
{
    int tkey = GetKey(layer, ladder);
    auto item = htable.find(tkey);
    return item != htable.end() ? item->second->size() : -1;
}

float HitTemporalIndexes::GetMinTime()
{
    float min_time { MAXTIME };
    for (auto item : htable)
    {
        min_time = min(min_time, item.second->top()->getTime());
    }
    return min_time;
}

float HitTemporalIndexes::GetMinTime(int layer, int ladder)
{
    int tkey = GetKey(layer, ladder);
    auto item = htable.find(tkey);
    if (item != htable.end()) return item->second->top()->getTime();
    return MAXTIME;
}

int HitTemporalIndexes::GetKey(int layer, int ladder)
{
    // TODO use cellID0 and/or cellID1
    return layer * 1000 + ladder;
}

float HitTemporalIndexes::MAXTIME { std::numeric_limits<float>::max() };
