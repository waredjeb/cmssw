// Authors: Marco Rovere - marco.rovere@cern.ch, Felice Pantaleo - felice.pantaleo@cern.ch
// Date: 09/2020

#ifndef RecoHGCal_TICL_ClusterFilterByAlgoAndSizeAndLayerRange_h
#define RecoHGCal_TICL_ClusterFilterByAlgoAndSizeAndLayerRange_h

#include "ClusterFilterBase.h"

#include <cassert>
#include <memory>
#include <utility>

// Filter clusters that belong to a specific algorithm
namespace ticl {
  class ClusterFilterByAlgoAndSizeAndLayerRange final : public ClusterFilterBase {
  public:
    ClusterFilterByAlgoAndSizeAndLayerRange(const edm::ParameterSet& ps)
        : ClusterFilterBase(ps),
          algo_number_(ps.getParameter<std::vector<int>>("algo_number")),
          min_cluster_size_(ps.getParameter<int>("min_cluster_size")),
          max_cluster_size_(ps.getParameter<int>("max_cluster_size")),
          min_layerId_(ps.getParameter<int>("min_layerId")),
          max_layerId_(ps.getParameter<int>("max_layerId")) {}
    ~ClusterFilterByAlgoAndSizeAndLayerRange() override {}

    void filter(const reco::CaloClusterSoAConstView& layerClusters,
                std::vector<float>& layerClustersMask,
                hgcal::RecHitTools& rhtools) const override {
      const int numberOfClusters = layerClusters.position().metadata().size();
      for (int i = 0; i < numberOfClusters; i++) {
        const DetId seedId = layerClusters.indexes()[i].seedID();
        assert(seedId.rawId() != 0);
        const auto layerId = rhtools.getLayerWithOffset(seedId);
        const auto cells = static_cast<unsigned int>(layerClusters.position()[i].cells());
        if (find(algo_number_.begin(), algo_number_.end(), layerClusters.indexes()[i].algoID()) ==
                algo_number_.end() or
            layerId > max_layerId_ or layerId < min_layerId_ or cells > max_cluster_size_ or
            (cells < min_cluster_size_ and rhtools.isSilicon(seedId))) {
          layerClustersMask[i] = 0.;
        }
      }
    }

  private:
    std::vector<int> algo_number_;
    unsigned int min_cluster_size_;
    unsigned int max_cluster_size_;
    unsigned int min_layerId_;
    unsigned int max_layerId_;
  };
}  // namespace ticl

#endif
