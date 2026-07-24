#ifndef RecoLocalCalo_HGCalRecProducers_plugins_alpaka_HGCalCLUEsteringAlgoWrapper_h
#define RecoLocalCalo_HGCalRecProducers_plugins_alpaka_HGCalCLUEsteringAlgoWrapper_h

#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // Device wrapper that runs the external, header-only CLUEstering library on
  // device to perform HGCal 2D per-layer layer-clustering. It is a drop-in
  // replacement for HGCalLayerClustersAlgoWrapper: same input/output SoA views,
  // same run() signature. Only clusterIndex(), isSeed() and
  // numberOfClustersScalar() of the output SoA are filled; delta/rho/
  // nearestHigher are left default-initialized (they are unused downstream).
  class HGCalCLUEsteringAlgoWrapper {
  public:
    void run(Queue& queue,
             const unsigned int size,
             const float dc,
             const float kappa,
             const float outlierDeltaFactor,
             const HGCalSoARecHitsDeviceCollection::ConstView inputs,
             HGCalSoARecHitsExtraDeviceCollection::View outputs) const;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoLocalCalo_HGCalRecProducers_plugins_alpaka_HGCalCLUEsteringAlgoWrapper_h
