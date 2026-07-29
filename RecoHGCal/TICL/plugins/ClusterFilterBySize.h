// Author: Marco Rovere - marco.rovere@cern.ch
// Date: 11/2018

#ifndef RecoHGCal_TICL_ClusterFilterBySize_h
#define RecoHGCal_TICL_ClusterFilterBySize_h

#include "ClusterFilterBase.h"

#include <memory>
#include <utility>

// Filter clusters that belong to a specific algorithm
namespace ticl {
  class ClusterFilterBySize final : public ClusterFilterBase {
  public:
    ClusterFilterBySize(const edm::ParameterSet& ps)
        : ClusterFilterBase(ps), max_cluster_size_(ps.getParameter<int>("max_cluster_size")) {}
    ~ClusterFilterBySize() override {}

    void filter(const reco::CaloClusterSoAConstView& layerClusters,
                std::vector<float>& layerClustersMask,
                hgcal::RecHitTools& rhtools) const override {
      const int numberOfClusters = layerClusters.position().metadata().size();
      for (int i = 0; i < numberOfClusters; i++) {
        if (static_cast<unsigned int>(layerClusters.position()[i].cells()) > max_cluster_size_) {
          layerClustersMask[i] = 0.;
        }
      }
    }

  private:
    unsigned int max_cluster_size_;
  };
}  // namespace ticl

#endif
