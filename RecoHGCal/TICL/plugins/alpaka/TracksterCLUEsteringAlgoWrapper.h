#ifndef RecoHGCal_TICL_plugins_alpaka_TracksterCLUEsteringAlgoWrapper_h
#define RecoHGCal_TICL_plugins_alpaka_TracksterCLUEsteringAlgoWrapper_h

#include <cstdint>
#include <span>
#include <vector>

#include "DataFormats/CaloRecHit/interface/alpaka/CaloClusterDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class TracksterCLUEsteringAlgoWrapper {
  public:
    // dc, ds, do in units sigmaT transversely and layerScale longitudinally
    struct Parameters {
      float dc;
      float rhoc;
      float outlierDistance;
      float seedingDistance;

      // Angular transverse scale [rad] per input collection (same order as `inputs` in run());
      // the density disc around a point has radius dc*sigmaT of its collection. Coordinates are
      // stored in units of the largest sigmaT and the metric rescales per pair, so showers
      // crossing a sub-detector boundary stay contiguous.
      std::vector<float> sigmaT;
      float layerScale;  // longitudinal scale in layers; the density cylinder has half-height dc*layerScale

      // Eta-dependent seed threshold: rhoc(point) = rhoc * (rhocPivotRadius / r)^rhocEtaExponent,
      // with r the point's gnomonic radius sqrt(x^2+y^2)/|z| (small r = high eta). 0 = flat rhoc.
      float rhocEtaExponent;
      float rhocPivotRadius;
    };

        void run(Queue& queue,
             std::span<const ::reco::CaloClusterSoAConstView> inputs,
             std::span<const uint32_t> offsets,
             const float* mask,
             uint32_t nTotal,
             uint32_t nSurviving,
             Parameters const& parameters,
             int32_t* assignment) const;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoHGCal_TICL_plugins_alpaka_TracksterCLUEsteringAlgoWrapper_h
