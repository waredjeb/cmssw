#ifndef RecoHGCal_TICL_interface_alpaka_CLUE3DParamsSoA_h
#define RecoHGCal_TICL_interface_alpaka_CLUE3DParamsSoA_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace ticl {

  // SoA for CLUE3D algorithm parameters
  // Supports up to 3 algorithm IDs (EM, HAD, SCINT typically)
  constexpr int kMaxAlgoIds = 3;

  GENERATE_SOA_LAYOUT(CLUE3DParamsSoALayout,
                      // Density calculation parameters (arrays indexed by algoId)
                      SOA_COLUMN(double, criticalDensity),
                      SOA_COLUMN(double, criticalSelfDensity),
                      SOA_COLUMN(int, densitySiblingLayers),
                      SOA_COLUMN(double, densityEtaPhiDistanceSqr),
                      SOA_COLUMN(double, densityXYDistanceSqr),
                      SOA_COLUMN(double, kernelDensityFactor),

                      // Seed finding parameters (arrays indexed by algoId)
                      SOA_COLUMN(double, criticalEtaPhiDistance),
                      SOA_COLUMN(double, criticalXYDistance),
                      SOA_COLUMN(int, criticalZDistanceLyr),
                      SOA_COLUMN(double, outlierMultiplier),
                      SOA_COLUMN(int, minNumLayerCluster),

                      // Boolean flags
                      SOA_SCALAR(bool, densityOnSameLayer),
                      SOA_SCALAR(bool, nearestHigherOnSameLayer),
                      SOA_SCALAR(bool, useAbsoluteProjectiveScale),
                      SOA_SCALAR(bool, useClusterDimensionXY),
                      SOA_SCALAR(bool, rescaleDensityByZ),

                      // Geometry info
                      SOA_SCALAR(int, lastLayerPerSide),  // number of layers on one side
                      SOA_SCALAR(int, nEtaBins),
                      SOA_SCALAR(int, nPhiBins),
                      SOA_SCALAR(float, etaMin),
                      SOA_SCALAR(float, etaMax),
                      SOA_SCALAR(float, phiMin),
                      SOA_SCALAR(float, phiMax)
  )

  using CLUE3DParamsSoA = CLUE3DParamsSoALayout<>;

}  // namespace ticl

#endif  // RecoHGCal_TICL_interface_alpaka_CLUE3DParamsSoA_h
