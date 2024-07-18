#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "RecoHGCal/TICL/interface/TracksterInferenceByDNN.h"
#include "RecoHGCal/TICL/interface/TracksterInferenceAlgoFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoHGCal/TICL/interface/PatternRecognitionAlgoBase.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "TrackstersPCA.h"

namespace ticl {
  using namespace cms::Ort;

  TracksterInferenceByDNN::TracksterInferenceByDNN(const edm::ParameterSet& conf)
    : TracksterInferenceAlgoBase(conf),
      eidInputName_(conf.getParameter<std::string>("eid_input_name")),
      eidOutputNameEnergy_(conf.getParameter<std::string>("eid_output_name_energy")),
      eidOutputNameId_(conf.getParameter<std::string>("eid_output_name_id")),
      eidMinClusterEnergy_(conf.getParameter<double>("eid_min_cluster_energy")),
      eidNLayers_(conf.getParameter<int>("eid_n_layers")),
      eidNClusters_(conf.getParameter<int>("eid_n_clusters"))
  {
    std::string id_modelPath = conf.getParameter<std::string>("RecoHGCal/TICL/data/onnx_models/id_v0.onnx");
    std::string en_modelPath = conf.getParameter<std::string>("RecoHGCal/TICL/data/onnx_models/energy_v0.onnx");

    static std::unique_ptr<cms::Ort::ONNXRuntime> onnxPIDRuntimeInstance = std::make_unique<cms::Ort::ONNXRuntime>(id_modelPath.c_str());
    onnxPIDSession_ = onnxPIDRuntimeInstance.get();

    static std::unique_ptr<cms::Ort::ONNXRuntime> onnxEnergyRuntimeInstance = std::make_unique<cms::Ort::ONNXRuntime>(en_modelPath.c_str());
    onnxEnergySession_ = onnxEnergyRuntimeInstance.get();

  }

  void TracksterInferenceByDNN::inputData(const std::vector<reco::CaloCluster>& layerClusters, std::vector<Trackster>& tracksters) {
    tracksterIndices.clear();
    for (int i = 0; i < static_cast<int>(tracksters.size()); i++) {
      float sumClusterEnergy = 0.;
      for (const unsigned int& vertex : tracksters[i].vertices()) {
        sumClusterEnergy += static_cast<float>(layerClusters[vertex].energy());
        if (sumClusterEnergy >= eidMinClusterEnergy_) {
          tracksters[i].setRegressedEnergy(0.f);
          tracksters[i].zeroProbabilities();
          tracksterIndices.push_back(i);
          break;
        }
      }
    }

    batchSize = static_cast<int>(tracksterIndices.size());
    if (batchSize == 0) return;

    std::vector<int64_t> inputShape = {batchSize, eidNLayers_, eidNClusters_, eidNFeatures_};
    input_shapes = {inputShape};

    input_Data.clear();
    input_Data.emplace_back(batchSize * eidNLayers_ * eidNClusters_ * eidNFeatures_, 0);

    inputNames = {"input:0"};
    outNames = {"output/regressed_energy:0", "output/id_probabilities:0"};

    for (int i = 0; i < batchSize; i++) {
      const Trackster& trackster = tracksters[tracksterIndices[i]];

      std::vector<int> clusterIndices(trackster.vertices().size());
      for (int k = 0; k < static_cast<int>(trackster.vertices().size()); k++) {
        clusterIndices[k] = k;
      }

      std::sort(clusterIndices.begin(), clusterIndices.end(), [&layerClusters, &trackster](const int& a, const int& b) {
        return layerClusters[trackster.vertices(a)].energy() > layerClusters[trackster.vertices(b)].energy();
      });

      std::vector<int> seenClusters(eidNLayers_, 0);

      for (const int& k : clusterIndices) {
        const reco::CaloCluster& cluster = layerClusters[trackster.vertices(k)];
        int j = rhtools_.getLayerWithOffset(cluster.hitsAndFractions()[0].first) - 1;
        if (j < eidNLayers_ && seenClusters[j] < eidNClusters_) {
          int index = (i * eidNLayers_ + j) * eidNClusters_ + seenClusters[j] * eidNFeatures_;
          input_Data[0][index] = static_cast<float>(cluster.energy() / static_cast<float>(trackster.vertex_multiplicity(k)));
          input_Data[0][index + 1] = static_cast<float>(std::abs(cluster.eta()));
          input_Data[0][index + 2] = static_cast<float>(cluster.phi());
          seenClusters[j]++;
        }
      }
    }
  }

  void TracksterInferenceByDNN::runInference(std::vector<Trackster>& tracksters) {

    if (!eidOutputNameEnergy_.empty()) {
      std::vector<float> energyOutputTensor = onnxEnergySession_->run(inputNames, input_Data, input_shapes, outNames, batchSize)[0];
      for (int i = 0; i < batchSize; i++) {
        float energy = energyOutputTensor[i];
        tracksters[tracksterIndices[i]].setRegressedEnergy(energy);
      }
    }

    if (!eidOutputNameId_.empty()) {
      std::vector<float> pidOutputTensor = onnxPIDSession_->run(inputNames, input_Data, input_shapes, outNames, batchSize)[0];
      float* probs = pidOutputTensor.data();
      for (int i = 0; i < batchSize; i++) {
        tracksters[tracksterIndices[i]].setProbabilities(probs);
        probs += tracksters[tracksterIndices[i]].id_probabilities().size();
        const auto& arr = tracksters[tracksterIndices[i]].id_probabilities();
      }
    }
  }
}

// Define this as a plug-in
DEFINE_EDM_PLUGIN(TracksterInferenceAlgoFactory, ticl::TracksterInferenceByDNN, "TracksterInferenceByDNN");
