// Alpaka port of PatternRecognitionbyCLUE3D
// Original author: Marco Rovere - marco.rovere@cern.ch

#ifndef RecoHGCal_TICL_PatternRecognitionByCLUE3D_Alpaka_H
#define RecoHGCal_TICL_PatternRecognitionByCLUE3D_Alpaka_H

#include "DataFormats/HGCalReco/interface/HGCalSoAClusters.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "RecoHGCal/TICL/interface/alpaka/PatternRecognitionAlgoBase.h"
#include <vector>
#include <memory>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class PatternRecognitionByCLUE3D final : public PatternRecognitionAlgoBase {
  public:
    PatternRecognitionByCLUE3D(const edm::ParameterSet& config);
    ~PatternRecognitionByCLUE3D() override = default;

    void makeTracksters(Queue& queue,
                        const HGCalSoAClustersDeviceCollection& layerClusters,
                        const TilesConstViewArray& tiles,
                        std::vector<ticl::Trackster>& tracksters) override;

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc);

  private:
    // Algorithm parameters from configuration
    const std::vector<double> criticalDensity_;
    const std::vector<double> criticalSelfDensity_;
    const std::vector<int> densitySiblingLayers_;
    const std::vector<double> densityEtaPhiDistanceSqr_;
    const std::vector<double> densityXYDistanceSqr_;
    const std::vector<double> kernelDensityFactor_;
    const bool densityOnSameLayer_;
    const bool nearestHigherOnSameLayer_;
    const bool useAbsoluteProjectiveScale_;
    const bool useClusterDimensionXY_;
    const bool rescaleDensityByZ_;
    const std::vector<double> criticalEtaPhiDistance_;
    const std::vector<double> criticalXYDistance_;
    const std::vector<int> criticalZDistanceLyr_;
    const std::vector<double> outlierMultiplier_;
    const std::vector<int> minNumLayerCluster_;
    const bool doPidCut_;
    const float cutHadProb_;
    const bool computeLocalTime_;
    const bool usePCACleaning_;

    // Device buffers for intermediate results
    // TODO: Define device buffers for:
    // - Local density calculation
    // - Nearest-higher distance calculation
    // - Seed finding
    // - Cluster assignment

    // Helper methods (to be implemented with kernels)
    // TODO: Implement these as Alpaka kernels
    // void calculateLocalDensity(Queue& queue, ...);
    // void calculateDistanceToHigher(Queue& queue, ...);
    // int findAndAssignTracksters(Queue& queue, ...);
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
