#pragma once

#include <limits>
#include <string>
#include <tuple>
#include <unordered_set>

#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"

#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

struct TrackerHitTimeWindowFilter final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::TrackerHitPlaneCollection,
                                            edm4hep::TrackerHitSimTrackerHitLinkCollection,
                                            edm4hep::TrackerHitSimTrackerHitLinkCollection>(
          const edm4hep::TrackerHitPlaneCollection&,
          const edm4hep::TrackerHitSimTrackerHitLinkCollection&,
          const edm4hep::TrackerHitSimTrackerHitLinkCollection&)>
{
  TrackerHitTimeWindowFilter(const std::string& name, ISvcLocator* svcLoc);

  std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection,
             edm4hep::TrackerHitSimTrackerHitLinkCollection>
  operator()(const edm4hep::TrackerHitPlaneCollection& hits,
             const edm4hep::TrackerHitSimTrackerHitLinkCollection& relations,
             const edm4hep::TrackerHitSimTrackerHitLinkCollection& rawHitsRelations) const override;

private:
  Gaudi::Property<float> m_timeWindowMin{this, "TimeWindowMin", -std::numeric_limits<float>::max(),
                                          "Minimum accepted hit time [ns] (inclusive)"};
  Gaudi::Property<float> m_timeWindowMax{this, "TimeWindowMax", std::numeric_limits<float>::max(),
                                          "Maximum accepted hit time [ns] (inclusive)"};
};
