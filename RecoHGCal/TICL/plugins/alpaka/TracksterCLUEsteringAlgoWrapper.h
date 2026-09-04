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
    struct Parameters {
      // dc, outlierDistance, seedingDistance in units of sigmaT
      float dc;
      float rhoc;
      float outlierDistance;
      float seedingDistance;

      // transversal angular width, one value per CE-E, Si CE-H and Sci CE-H
      std::vector<float> sigmaT;

      float layerScale;  // longitudinal scale in layers

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
