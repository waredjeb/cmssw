#ifndef RecoLocalCalo_HGCalRecProducers_LayerClustersAndAssociations_h
#define RecoLocalCalo_HGCalRecProducers_LayerClustersAndAssociations_h

#include "DataFormats/CaloRecHit/interface/CaloClusterHostCollection.h"
#include "DataFormats/TICL/interface/HitsAndFractionsHost.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"

#include <memory>

namespace ticl {
  // Return the portable SoA of per-cluster information, plus variable length
  // hit lists keyed by cluster index in CSR format.
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
