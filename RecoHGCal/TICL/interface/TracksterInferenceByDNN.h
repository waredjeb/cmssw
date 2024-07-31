#ifndef RecoHGCal_TICL_TracksterInferenceByDNN_H__
#define RecoHGCal_TICL_TracksterInferenceByDNN_H__

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "RecoHGCal/TICL/interface/TracksterInferenceAlgoBase.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

namespace ticl {
  using namespace cms::Ort;

  class TracksterInferenceByDNN : public TracksterInferenceAlgoBase {
  public:
    explicit TracksterInferenceByDNN(const edm::ParameterSet& conf);
    void inputData(const std::vector<reco::CaloCluster>& layerClusters, std::vector<Trackster>& tracksters) override;
    void runInference(std::vector<Trackster>& tracksters, const std::string& mode, const std::string& operation) override;

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc);

  private:
    const cms::Ort::ONNXRuntime* onnxPIDSessionCLU3D_;
    const cms::Ort::ONNXRuntime* onnxEnergySessionCLU3D_;
    const cms::Ort::ONNXRuntime* onnxPIDSessionLinking_;
    const cms::Ort::ONNXRuntime* onnxEnergySessionLinking_;

    std::string id_modelPath_CLU3D;
    std::string en_modelPath_CLU3D;
    std::string id_modelPath_Linking;
    std::string en_modelPath_Linking;

    //const std::string eidInputName_;
    const std::string eidOutputNameEnergy_;
    const std::string eidOutputNameId_;
    const float eidMinClusterEnergy_;
    const int eidNLayers_;
    const int eidNClusters_;
    static const int eidNFeatures_ = 3;

    hgcal::RecHitTools rhtools_;
    std::vector<std::vector<int64_t>> input_shapes;
    std::vector<int> tracksterIndices;
    std::vector<std::vector<float>> input_Data;
    int batchSize;
  };
}  // namespace ticl

#endif  // RecoHGCal_TICL_TracksterInferenceByDNN_H__
