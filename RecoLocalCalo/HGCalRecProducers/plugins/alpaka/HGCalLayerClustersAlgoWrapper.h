#ifndef RecoLocalCalo_HGCalRecProducers_plugins_alpaka_HGCalLayerClustersAlgoWrapper_h
#define RecoLocalCalo_HGCalRecProducers_plugins_alpaka_HGCalLayerClustersAlgoWrapper_h

#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include <cstdint>
#include <span>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class HGCalLayerClustersAlgoWrapper {
  public:
    // batchItemSizes has number of hits for each non empty layer, allowing
    // independent clustering of layers
    void run(Queue& queue,
             const unsigned int size,
             const float dc,
             const float kappa,
             const float outlierDeltaFactor,
             std::span<const uint32_t> batchItemSizes,
             const HGCalSoARecHitsDeviceCollection::ConstView inputs,
             HGCalSoARecHitsExtraDeviceCollection::View outputs) const;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoLocalCalo_HGCalRecProducers_plugins_alpaka_HGCalLayerClustersAlgoWrapper_h
