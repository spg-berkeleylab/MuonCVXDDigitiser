#include "TrackerHitTimeWindowFilter.h"

#include <cmath>

#include "TMath.h"

namespace {
// Same propagation-time correction DDPlanarDigi applies before its own
// UseTimeWindow check (hitT -= r/c, where r is the full 3D distance from
// the origin) - see DDPlanarDigi.cpp. MuonCVXDDigitiser does not apply
// this correction itself, so we do it here, before filtering, to match
// DDPlanarDigi's ordering (correct, then window-check) instead of
// filtering on the raw digitized time.
double correctedTime(double time, const edm4hep::Vector3d& pos) {
  const double r = std::sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);
  const double dt = r / (TMath::C() / 1e6);
  return time - dt;
}
}  // namespace

TrackerHitTimeWindowFilter::TrackerHitTimeWindowFilter(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(
          name, svcLoc,
          {KeyValues("InputHits", {"VXDBarrelHits"}),
           KeyValues("InputRelations", {"VXDBarrelHitsRelations"}),
           KeyValues("InputRawHitsRelations", {"VXDBarrelRawHitsRelations"})},
          {KeyValues("OutputHits", {"VXDBarrelHits_TimeFiltered"}),
           KeyValues("OutputRelations", {"VXDBarrelHitsRelations_TimeFiltered"}),
           KeyValues("OutputRawHitsRelations", {"VXDBarrelRawHitsRelations_TimeFiltered"})})
{
}

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection,
           edm4hep::TrackerHitSimTrackerHitLinkCollection>
TrackerHitTimeWindowFilter::operator()(const edm4hep::TrackerHitPlaneCollection& hits,
                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& relations,
                                        const edm4hep::TrackerHitSimTrackerHitLinkCollection& rawHitsRelations) const
{
  edm4hep::TrackerHitPlaneCollection outHits;
  outHits.setSubsetCollection();

  // Track which hits (by their index in `hits`) survived the window, so
  // the relation collections below can be filtered to match.
  std::unordered_set<int> keptIndices;

  for (const auto& hit : hits) {
    const double t = correctedTime(hit.getTime(), hit.getPosition());
    if (t >= m_timeWindowMin.value() && t <= m_timeWindowMax.value()) {
      outHits.push_back(hit);
      keptIndices.insert(hit.getObjectID().index);
    }
  }

  edm4hep::TrackerHitSimTrackerHitLinkCollection outRelations;
  outRelations.setSubsetCollection();
  for (const auto& rel : relations) {
    const auto from = rel.getFrom();
    if (from.isAvailable() && keptIndices.count(from.getObjectID().index)) {
      outRelations.push_back(rel);
    }
  }

  edm4hep::TrackerHitSimTrackerHitLinkCollection outRawHitsRelations;
  outRawHitsRelations.setSubsetCollection();
  for (const auto& rel : rawHitsRelations) {
    const auto from = rel.getFrom();
    if (from.isAvailable() && keptIndices.count(from.getObjectID().index)) {
      outRawHitsRelations.push_back(rel);
    }
  }

  debug() << "TrackerHitTimeWindowFilter: kept " << outHits.size() << "/" << hits.size() << " hits, "
          << outRelations.size() << "/" << relations.size() << " relations, " << outRawHitsRelations.size() << "/"
          << rawHitsRelations.size() << " raw-hits relations, in window [" << m_timeWindowMin.value() << ", "
          << m_timeWindowMax.value() << "] ns" << endmsg;

  return std::make_tuple(std::move(outHits), std::move(outRelations), std::move(outRawHitsRelations));
}

DECLARE_COMPONENT(TrackerHitTimeWindowFilter)
