// Author: Wahid Redjeb - wahid.redjeb@cern.ch
// Date: 01/2025

#ifndef RecoHGCal_TICL_TracksterFilterByPDGID_h
#define RecoHGCal_TICL_TracksterFilterByPDGID_h

#include "RecoHGCal/TICL/plugins/TracksterFilterBase.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"

#include <memory>
#include <utility>

// Filter tracksters based on their hadronic/EM nature
namespace ticl {
  class TracksterFilterByPDGID final : public TracksterFilterBase {
  public:
    TracksterFilterByPDGID(const edm::ParameterSet& ps)
        : TracksterFilterBase(ps), keep_hadronic_(ps.getParameter<bool>("keepHadronic")) {}
    ~TracksterFilterByPDGID() override = default;

    void filter(const std::vector<ticl::Trackster>& tracksters,
                const std::vector<reco::CaloCluster>& layerClusters,
                std::vector<float>& trackstersMask,
                hgcal::RecHitTools& rhtools) const override {
      for (size_t i = 0; i < tracksters.size(); ++i) {
        if (trackstersMask[i] == 0.f)
          continue;  // Already masked

        const bool is_hadronic = tracksters[i].isHadronic();

        // Mask tracksters that don't match criterion
        if (keep_hadronic_ && !is_hadronic) {
          trackstersMask[i] = 0.f;
        } else if (!keep_hadronic_ && is_hadronic) {
          trackstersMask[i] = 0.f;
        }
      }
    }

  private:
    bool keep_hadronic_;
  };
}  // namespace ticl

#endif
