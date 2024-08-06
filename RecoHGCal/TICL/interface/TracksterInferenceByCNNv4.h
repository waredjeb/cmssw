#ifndef RecoHGCal_TICL_TracksterInferenceByCNNv4_H__
#define RecoHGCal_TICL_TracksterInferenceByCNNv4_H__

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "RecoHGCal/TICL/interface/TracksterInferenceAlgoBase.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

namespace ticl {
  using namespace cms::Ort;

  class TracksterInferenceByCNNv4 : public TracksterInferenceAlgoBase {
  public:
    explicit TracksterInferenceByCNNv4(const edm::ParameterSet& conf);
    void inputData(const std::vector<reco::CaloCluster>& layerClusters, std::vector<Trackster>& tracksters) override;
    void runInference(std::vector<Trackster>& tracksters) override;

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc);
  private:
    const cms::Ort::ONNXRuntime* onnxSession_;

    const std::string modelPath_;
    const std::string eidOutputNameEnergy_;
    const std::string eidOutputNameId_;
    const float eidMinClusterEnergy_;
    const int eidNLayers_;
    const int eidNClusters_;
    static constexpr int eidNFeatures_ = 3;
    int doPID_;
    int doRegression_;

    hgcal::RecHitTools rhtools_;
    std::vector<std::vector<int64_t>> input_shapes;
    std::vector<int> tracksterIndices;
    std::vector<std::vector<float>> input_Data;
    int batchSize;
  };
}  // namespace ticl

#endif  // RecoHGCal_TICL_TracksterInferenceByDNN_H__
