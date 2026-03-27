
#pragma once

#include "DataFormats/TICL/interface/CaloClusterHostCollection.h"
#include "DataFormats/TICL/interface/HitsAndFractionsHost.h"
#include <memory>

namespace ticl {

  struct LayerClustersAndAssociations {
    std::unique_ptr<reco::CaloClusterHostCollection> layer_clusters;
    // TESTING: Removed to isolate crash
    // std::unique_ptr<HitsAndFractionsHost> hits_and_fractions;

    LayerClustersAndAssociations(int number_of_clusters, int total_rechits)
        : layer_clusters{std::make_unique<reco::CaloClusterHostCollection>(
              cms::alpakatools::host(), number_of_clusters, number_of_clusters, number_of_clusters, number_of_clusters)} {
          // Don't construct hits_and_fractions for testing
    }
  };

}  // namespace ticl
