#ifndef RecoHGCal_TICL_plugins_alpaka_CLUE3DKernel_h
#define RecoHGCal_TICL_plugins_alpaka_CLUE3DKernel_h

#include "DataFormats/HGCalReco/interface/alpaka/CLUE3DStateDeviceCollection.h"
#include "RecoHGCal/TICL/plugins/alpaka/HGCalTiles.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoHGCal/TICL/interface/alpaka/CLUE3DParamsSoA.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace ticl {

    class CLUE3DKernel {
    public:
      CLUE3DKernel() = default;

      // Calculate local density for all clusters
      void calculateLocalDensity(Queue& queue,
                                  const ::ticl::CLUE3DParamsSoA::ConstView params,
                                  const HGCalTilesView& tiles,
                                  const float* layersPosZ,  // array of Z positions per layer
                                  CLUE3DStateDeviceCollection& state,
                                  int nClusters);

      // Calculate distance to nearest higher density cluster
      void calculateDistanceToHigher(Queue& queue,
                                      const ::ticl::CLUE3DParamsSoA::ConstView params,
                                      const HGCalTilesView& tiles,
                                      CLUE3DStateDeviceCollection& state,
                                      int nClusters);

      // Find seeds and classify outliers
      void findAndAssignSeeds(Queue& queue,
                              const ::ticl::CLUE3DParamsSoA::ConstView params,
                              CLUE3DStateDeviceCollection& state,
                              int nClusters,
                              int* nSeeds);  // output: number of seeds found
    };

  }  // namespace ticl
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoHGCal_TICL_plugins_alpaka_CLUE3DKernel_h
