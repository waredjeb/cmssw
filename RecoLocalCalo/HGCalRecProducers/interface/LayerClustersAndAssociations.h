#ifndef RecoLocalCalo_HGCalRecProducers_LayerClustersAndAssociations_h
#define RecoLocalCalo_HGCalRecProducers_LayerClustersAndAssociations_h

#include "DataFormats/CaloRecHit/interface/CaloClusterHostCollection.h"
#include "DataFormats/TICL/interface/HitsAndFractionsHost.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"

#include <memory>

namespace ticl {

  // Return bundle of the layer-cluster algorithms: the per-cluster scalars as a
  // portable SoA, plus the variable-length hit lists as a separate compressed
  // (CSR) association map keyed by cluster index. The two are index-consistent:
  // hits_and_fractions[i] holds the hits of cluster i of layer_clusters.
  struct LayerClustersAndAssociations {
    std::unique_ptr<reco::CaloClusterHostCollection> layer_clusters;
    std::unique_ptr<HitsAndFractionsHost> hits_and_fractions;

    LayerClustersAndAssociations(int number_of_clusters, int total_rechits)
        : layer_clusters{std::make_unique<reco::CaloClusterHostCollection>(
              cms::alpakatools::host(), number_of_clusters, number_of_clusters, number_of_clusters, number_of_clusters)},
          hits_and_fractions{
              std::make_unique<HitsAndFractionsHost>(cms::alpakatools::host(), total_rechits, number_of_clusters)} {}
  };

}  // namespace ticl

#endif
