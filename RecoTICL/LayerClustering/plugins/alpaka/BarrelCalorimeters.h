#ifndef RecoTICL_LayerClustering_plugins_alpaka_BarrelCalorimeters_h
#define RecoTICL_LayerClustering_plugins_alpaka_BarrelCalorimeters_h

#include <cstdint>

#include "DataFormats/ParticleFlowReco/interface/alpaka/PFRecHitDeviceCollection.h"
#include "CondFormats/DataRecord/interface/EcalPFRecHitThresholdsRcd.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoParticleFlow/PFRecHitProducer/interface/PFRecHitTopologyRecord.h"
#include "RecoParticleFlow/PFRecHitProducer/interface/alpaka/PFRecHitParamsDeviceCollection.h"
#include "RecoParticleFlow/PFRecHitProducer/interface/alpaka/PFRecHitTopologyDeviceCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::barrelclue {

  struct ECAL {
    static constexpr uint32_t kNumberOfLayers = 1;

    ALPAKA_FN_HOST_ACC static uint32_t layer(reco::PFRecHitDeviceCollection::ConstView hits, uint32_t i) { return 0u; }

    using NoiseType = reco::PFRecHitECALParamsDeviceCollection;
    using NoiseRecordType = EcalPFRecHitThresholdsRcd;

    ALPAKA_FN_HOST_ACC static float sigmaNoise(NoiseType::ConstView noise, uint32_t denseId) {
      return noise[denseId].energyThresholds();
    }
  };

  struct HCAL {
    static constexpr uint32_t kNumberOfLayers = 8;

    ALPAKA_FN_HOST_ACC static uint32_t layer(reco::PFRecHitDeviceCollection::ConstView hits, uint32_t i) {
      return static_cast<uint32_t>(hits[i].depth());
    }

    using NoiseType = reco::PFRecHitHCALTopologyDeviceCollection;
    using NoiseRecordType = PFRecHitHCALTopologyRecord;

    ALPAKA_FN_HOST_ACC static float sigmaNoise(NoiseType::ConstView noise, uint32_t denseId) {
      // The thresholds are only filled when the ESProducer runs with usePFThresholdsFromDB = True
      return noise.cutsFromDB() ? noise.seedThreshold()[denseId] : 0.f;
    }
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::barrelclue

#endif  // RecoTICL_LayerClustering_plugins_alpaka_BarrelCalorimeters_h
