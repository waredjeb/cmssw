// Author: Wahid Redjeb - wahid.wahid@cern.ch
// Date: 06/2025

#ifndef RecoHGCal_TICL_TracksterFilterBySize_H__
#define RecoHGCal_TICL_TracksterFilterBySize_H__

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "TracksterFilterBase.h"

#include <memory>
#include <utility>

// Filter clusters that belong to a specific algorithm
namespace ticl {
  class TracksterFilterBySize final : public TracksterFilterBase {
  public:
    TracksterFilterBySize(const edm::ParameterSet& ps)
        : TracksterFilterBase(ps), max_cluster_size_(ps.getParameter<int>("max_cluster_size")) {}
    ~TracksterFilterBySize() override {}

    void filter(const std::vector<ticl::Trackster>& tracksters,
                const std::vector<reco::CaloCluster>& layerClusters,
                std::vector<float>& trackstersMask,
                hgcal::RecHitTools& rhtools) const override {
      for (size_t i = 0; i < tracksters.size(); i++) {
        if (tracksters[i].vertices().size() > max_cluster_size_) {
          trackstersMask[i] = 0.f;
        }
      }
    }

  private:
    unsigned int max_cluster_size_;
  };
}  // namespace ticl

#endif
