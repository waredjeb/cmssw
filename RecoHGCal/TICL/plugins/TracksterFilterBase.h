// Author: Wahid Redjeb - wahid.redjeb@cern.ch
// Date: 01/2025

#ifndef RecoHGCal_TICL_TracksterFilterBase_h
#define RecoHGCal_TICL_TracksterFilterBase_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include <vector>

namespace edm {
  class ParameterSet;
}
namespace reco {
  class CaloCluster;
}

namespace ticl {
  class TracksterFilterBase {
  public:
    explicit TracksterFilterBase(const edm::ParameterSet&) {}
    virtual ~TracksterFilterBase() = default;

    virtual void filter(const std::vector<ticl::Trackster>& tracksters,
                        const std::vector<reco::CaloCluster>& layerClusters,
                        std::vector<float>& trackstersMask,
                        hgcal::RecHitTools& rhtools) const = 0;
  };
}  // namespace ticl

#endif
