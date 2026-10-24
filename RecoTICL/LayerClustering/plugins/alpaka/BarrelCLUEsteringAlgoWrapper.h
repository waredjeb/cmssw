#ifndef RecoTICL_LayerClustering_plugins_alpaka_BarrelCLUEsteringAlgoWrapper_h
#define RecoTICL_LayerClustering_plugins_alpaka_BarrelCLUEsteringAlgoWrapper_h

#include <cstdint>
#include <span>

// PFRecHitDeviceCollection pulls in Eigen, which must precede any other SoA header
#include "DataFormats/ParticleFlowReco/interface/alpaka/PFRecHitDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

#include "BarrelCalorimeters.h"

namespace barrelclue {
  GENERATE_SOA_LAYOUT(BarrelCLUEsteringInputsLayout,
                      SOA_COLUMN(float, eta),
                      SOA_COLUMN(float, phi),
                      SOA_COLUMN(float, energy),
                      SOA_COLUMN(float, sigmaNoise),
                      SOA_COLUMN(uint32_t, detId),
                      SOA_COLUMN(uint32_t, hitIndex),
                      SOA_COLUMN(int32_t, clusterIndex))

  using BarrelCLUEsteringInputsSoA = BarrelCLUEsteringInputsLayout<>;
}  // namespace barrelclue

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using BarrelCLUEsteringInputsDeviceCollection = PortableCollection<::barrelclue::BarrelCLUEsteringInputsSoA>;

  template <typename CAL>
  class BarrelCLUEsteringAlgoWrapper {
  public:
    void countHitsPerLayer(Queue& queue,
                           const reco::PFRecHitDeviceCollection::ConstView inputs,
                           std::span<uint32_t> layerSizes) const;

    void run(Queue& queue,
             uint32_t size,
             float dc,
             float kappa,
             std::span<const uint32_t> layerSizes,
             const reco::PFRecHitDeviceCollection::ConstView inputs,
             const typename CAL::NoiseType::ConstView noise,
             HGCalSoARecHitsExtraDeviceCollection::View outputs) const;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoTICL_LayerClustering_plugins_alpaka_BarrelCLUEsteringAlgoWrapper_h
