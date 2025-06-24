// Author: Wahid Redjeb - wahid.redjeb@cern.ch 
// Date: 06/2025

#ifndef RecoHGCal_TICL_TracksterFilterByAlgo_H__
#define RecoHGCal_TICL_TracksterFilterByAlgo_H__

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "TracksterFilterBase.h"
#include <memory>
#include <utility>

// Filter clusters that belong to a specific algorithm
namespace ticl {
  class TracksterFilterByAlgo final : public TracksterFilterBase {
  public:
    TracksterFilterByAlgo(const edm::ParameterSet& ps)
        : TracksterFilterBase(ps), algo_number_(ps.getParameter<std::vector<int>>("algo_number")) {}
    ~TracksterFilterByAlgo() override {}

    void filter(const std::vector<ticl::Trackster>& tracksters,
                const std::vector<reco::CaloCluster>& layerClusters,
                std::vector<float>& tracksterMask,
                hgcal::RecHitTools& rhtools) const override {
      for (size_t i = 0; i < tracksters.size(); i++) {
        if (find(algo_number_.begin(), algo_number_.end(), layerClusters[i].algo()) == algo_number_.end()) {
          tracksterMask[i] = 0.;
        }
      }
    }

  private:
    std::vector<int> algo_number_;
  };
}  // namespace ticl

#endif
