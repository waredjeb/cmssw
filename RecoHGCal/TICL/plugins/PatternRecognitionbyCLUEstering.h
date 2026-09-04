#ifndef RecoHGCal_TICL_PatternRecognitionbyCLUEstering_h
#define RecoHGCal_TICL_PatternRecognitionbyCLUEstering_h

#include <cstdint>
#include <vector>

#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "RecoHGCal/TICL/interface/PatternRecognitionAlgoBase.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/TICLGeomTools.h"

namespace ticl {

  // Builds tracksters from the layer cluster to trackster assignment computed on device
  // This plugin groups the layer clusters and fill the trackster information, but clustering
  // has already happened.
  template <typename TILES>
  class PatternRecognitionbyCLUEstering final : public PatternRecognitionAlgoBaseT<TILES> {
  public:
    PatternRecognitionbyCLUEstering(const edm::ParameterSet& conf, edm::ConsumesCollector);
    ~PatternRecognitionbyCLUEstering() override = default;

    void makeTracksters(const typename PatternRecognitionAlgoBaseT<TILES>::Inputs& input,
                        std::vector<Trackster>& result,
                        std::unordered_map<int, std::vector<int>>& seedToTracksterAssociation) override;

    void filter(std::vector<Trackster>& output,
                const std::vector<Trackster>& inTracksters,
                const typename PatternRecognitionAlgoBaseT<TILES>::Inputs& input,
                std::unordered_map<int, std::vector<int>>& seedToTracksterAssociation) override;

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc);
    void setGeometry(ticlgeom::Tools const& rhtools) override;

    // Masks the layer clusters CLUE classified as outliers (noise), so that the Recovery iteration
    // does not turn them into one-cluster tracksters.
    void maskLayerClusters(const typename PatternRecognitionAlgoBaseT<TILES>::Inputs& input,
                           std::vector<float>& outputMask) const override;

  private:
    const edm::EDGetTokenT<std::vector<int32_t>> assignmentToken_;
    const bool maskOutliers_;
    const bool doPidCut_;
    const float cutHadProb_;
    const bool computeLocalTime_;
    const bool usePCACleaning_;
    const int minNumLayerCluster_;
  };

}  // namespace ticl

#endif  // RecoHGCal_TICL_PatternRecognitionbyCLUEstering_h
