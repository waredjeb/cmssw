#pragma once

#include "DataFormats/HGCalReco/interface/HGCalSoAClusters.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoHGCal/TICL/interface/alpaka/PatternRecognitionAlgoBase.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/HGCalReco/interface/Tiles.h"
#include "DataFormats/HGCalReco/interface/alpaka/TilesDevice.h"

#include <algorithm>
#include <array>
#include <vector>
#include <unordered_map>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class PatternRecognitionByCLUE3D final : public PatternRecognitionAlgoBase {
  private:
    // Algorithm parameters
    std::vector<double> criticalDensity_;
    std::vector<double> criticalSelfDensity_;
    std::vector<int> densitySiblingLayers_;
    std::vector<double> densityEtaPhiDistanceSqr_;
    std::vector<double> densityXYDistanceSqr_;
    std::vector<double> kernelDensityFactor_;
    bool densityOnSameLayer_;
    bool nearestHigherOnSameLayer_;
    bool useAbsoluteProjectiveScale_;
    bool useClusterDimensionXY_;
    bool rescaleDensityByZ_;
    std::vector<double> criticalEtaPhiDistance_;
    std::vector<double> criticalXYDistance_;
    std::vector<int> criticalZDistanceLyr_;
    std::vector<double> outlierMultiplier_;
    std::vector<int> minNumLayerCluster_;
    bool doPidCut_;
    float cutHadProb_;
    bool computeLocalTime_;
    bool usePCACleaning_;

  public:
    PatternRecognitionByCLUE3D(const edm::ParameterSet& config);
    ~PatternRecognitionByCLUE3D() override = default;

    void makeTracksters(Queue& queue,
                        const HGCalSoAClustersDeviceCollection& layerClusters,
                        std::vector<::ticl::Trackster>& tracksters,
                        std::array<ticl::TICLLayerTilesDevice::View, 96> tiles) override;

    static void fillPSetDescription(::edm::ParameterSetDescription& iDesc);
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
