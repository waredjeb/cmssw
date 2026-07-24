#ifndef RecoLocalCalo_HGCalRecProducers_LayerClusterAndAssociations_h
#define RecoLocalCalo_HGCalRecProducers_LayerClusterAndAssociations_h

#include "DataFormats/CaloRecHit/interface/CaloClusterHostCollection.h"
#include "DataFormats/TICL/interface/HitsAndFractionsHost.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include <memory>

namespace ticl {

  // Bundle produced by the layer-clustering algorithms: the per-cluster SoA
  // (reco::CaloClusterHostCollection) plus the transient hits-and-fractions
  // association map (layer-cluster index -> {DetId, fraction}). The map is used
  // within the job (e.g. to compute cluster timing) and to repopulate the legacy
  // reco::CaloCluster::hitsAndFractions in the SoA -> legacy converter; it is not
  // persisted to the event.
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
