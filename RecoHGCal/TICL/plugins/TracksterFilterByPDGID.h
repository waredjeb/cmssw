// Author: Wahid Redjeb - wahid.wahid@cern.ch
// Date: 06/2025

#ifndef RecoHGCal_TICL_TracksterFilterByPDGID_H__
#define RecoHGCal_TICL_TracksterFilterByPDGID_H__

#include <oneapi/tbb/parallel_pipeline.h>
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "TracksterFilterBase.h"

#include <memory>
#include <ranges>
#include <utility>

// Filter clusters that belong to a specific algorithm
namespace ticl {
  class TracksterFilterByPDGID final : public TracksterFilterBase {
  public:
    TracksterFilterByPDGID(const edm::ParameterSet& ps)
        : TracksterFilterBase(ps),
          filterEM_(ps.getParameter<bool>("filterEM")),
          threshold_(ps.getParameter<double>("threshold")) {};
    ~TracksterFilterByPDGID() override {}

    void filter(const std::vector<ticl::Trackster>& tracksters,
                const std::vector<reco::CaloCluster>& layerClusters,
                std::vector<float>& trackstersMask,
                hgcal::RecHitTools& rhtools) const override {
      auto isEM = [this](const Trackster& t) -> bool {
        auto const emProb = t.id_probability(ticl::Trackster::ParticleType::electron) +
                            t.id_probability(ticl::Trackster::ParticleType::photon);
        return emProb >= threshold_;
      };

      std::ranges::transform(
          tracksters, trackstersMask.begin(), [&](const Trackster& t) { return (filterEM_ && isEM(t)) ? 0.f : 1.f; });
    }

  private:
    bool filterEM_;
    double threshold_;
  };
}  // namespace ticl

#endif
