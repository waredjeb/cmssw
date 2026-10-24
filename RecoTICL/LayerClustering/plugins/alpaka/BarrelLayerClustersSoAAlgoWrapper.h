#ifndef RecoTICL_LayerClustering_plugins_alpaka_BarrelLayerClustersSoAAlgoWrapper_h
#define RecoTICL_LayerClustering_plugins_alpaka_BarrelLayerClustersSoAAlgoWrapper_h

// PFRecHitDeviceCollection pulls in Eigen, which must precede any other SoA header
#include "DataFormats/ParticleFlowReco/interface/alpaka/PFRecHitDeviceCollection.h"
#include "DataFormats/CaloRecHit/interface/alpaka/CaloClusterDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include "BarrelCalorimeters.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  template <typename CAL>
  class BarrelLayerClustersSoAAlgoWrapper {
  public:
    void run(Queue& queue,
             const unsigned int numberOfClusters,
             const reco::PFRecHitDeviceCollection::ConstView inputRecHits,
             const HGCalSoARecHitsExtraDeviceCollection::ConstView inputClusters,
             reco::CaloClusterDeviceCollection::View outputs) const;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoTICL_LayerClustering_plugins_alpaka_BarrelLayerClustersSoAAlgoWrapper_h
