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
      id_modelPath(conf.getParameter<edm::FileInPath>("onnxPIDModelPath").fullPath()),
      en_modelPath(conf.getParameter<edm::FileInPath>("onnxEnergyModelPath").fullPath()), 
      eidMinClusterEnergy_(conf.getParameter<double>("eid_min_cluster_energy")),
      eidNLayers_(conf.getParameter<int>("eid_n_layers")),
      eidNClusters_(conf.getParameter<int>("eid_n_clusters"))
  {
    //std::string id_modelPath = conf.getParameter<std::string>("RecoHGCal/TICL/data/onnx_models/id_v0.onnx");
    //std::string en_modelPath = conf.getParameter<std::string>("RecoHGCal/TICL/data/onnx_models/energy_v0.onnx");

    //std::string id_modelPath = conf.getParameter<std::string>("RecoHGCal/TICL/data/tf_models/id_v0.onnx");
    //std::string en_modelPath = conf.getParameter<std::string>("RecoHGCal/TICL/data/tf_models/energy_v0.onnx");

    //desc.add<edm::FileInPath>("onnxPIDModelPath", edm::FileInPath("RecoHGCal/TICL/data/tf_models/id_v0.onnx"))
    //->setComment("Path ONNX PID model");
    //desc.add<edm::FileInPath>("onnxEnergyModelPath", edm::FileInPath("RecoHGCal/TICL/data/tf_models/energy_v0.onnx"))
    //->setComment("Path ONNX Enegry model");

  
    static std::unique_ptr<cms::Ort::ONNXRuntime> onnxPIDRuntimeInstance = std::make_unique<cms::Ort::ONNXRuntime>(id_modelPath.c_str());
    onnxPIDSession_ = onnxPIDRuntimeInstance.get();

    static std::unique_ptr<cms::Ort::ONNXRuntime> onnxEnergyRuntimeInstance = std::make_unique<cms::Ort::ONNXRuntime>(en_modelPath.c_str());
    onnxEnergySession_ = onnxEnergyRuntimeInstance.get();

  }

  void TracksterInferenceByDNN::inputData(const std::vector<reco::CaloCluster>& layerClusters, std::vector<Trackster>& tracksters) {
    //tracksterIndices.clear();
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
    if (batchSize == 0) return;

    std::vector<std::string> inputNames = {"input:0"};
    std::vector<std::string> output_e = {"output/regressed_energy:0"};//, "output/id_probabilities:0"};
    std::vector<std::string> output_i = {"output/id_probabilities:0"};

    //std::cout << "[DEBUG] inputNames = " << inputNames << std::endl;
    //std::cout << "[DEBUG] inputNames = " << inputNames.size() << std::endl;
    //std::cout << "[DEBUG] inputData = " << input_Data.size() << std::endl;
    //std::cout << "[DEBUG] batchSize = " << batchSize << std::endl;
    //std::cout << "[DEBUG] eidNLayers_ = " << eidNLayers_ << std::endl;
    //std::cout << "[DEBUG] eidNClusters_ = " << eidNClusters_ << std::endl;
    //std::cout << "[DEBUG] eidNFeatures_ = " << eidNFeatures_ << std::endl;
    // Print inputNames
    //std::cout << "inputNames:" << std::endl;
    //for (const auto& name : inputNames) {
    //std::cout << name << std::endl;
    //}

    // Print output_i
    //std::cout << "output_i:" << std::endl;
    //for (const auto& name : output_i) {
    //std::cout << name << std::endl;
    //}
    
    // Print output_e
    //std::cout << "output_e:" << std::endl;
    //for (const auto& name : output_e) {
    //std::cout << name << std::endl;
    //}
    //std::cout << "Printing input_Data:" << std::endl;
    //for (size_t i = 0; i < input_Data.size(); ++i) {
    //std::cout << "Row " << i << ": ";
    //for (size_t j = 0; j < input_Data[i].size(); ++j) {
    //      std::cout << input_Data[i][j] << " ";
    //}
    //std::cout << std::endl;
    //}
    
    std::vector<float> energyOutputTensor = onnxEnergySession_->run(inputNames, input_Data, input_shapes, output_e, batchSize)[0];
    //std::cout<< "size of outputTensors = " << energyOutputTensor.size() << std::endl;

    if (!output_e.empty()) {
      for (int i = 0; i < batchSize; i++) {
        float energy = energyOutputTensor[i];
        tracksters[tracksterIndices[i]].setRegressedEnergy(energy);
      }
    }
    std::vector<float> pidOutputTensor = onnxPIDSession_->run(inputNames, input_Data, input_shapes, output_i, batchSize)[0];
    
    if (!output_i.empty()) {
      float* probs = pidOutputTensor.data();
      for (int i = 0; i < batchSize; i++) {
        tracksters[tracksterIndices[i]].setProbabilities(probs);
        probs += tracksters[tracksterIndices[i]].id_probabilities().size();
      }
    }
  }
  void TracksterInferenceByDNN::fillPSetDescription(edm::ParameterSetDescription& iDesc) {
    iDesc.add<int>("algo_verbosity", 0);
    iDesc.add<edm::FileInPath>("onnxPIDModelPath", edm::FileInPath("RecoHGCal/TICL/data/tf_models/id_v0.onnx"))
      ->setComment("Path ONNX PID model");
    iDesc.add<edm::FileInPath>("onnxEnergyModelPath", edm::FileInPath("RecoHGCal/TICL/data/tf_models/energy_v0.onnx"))
      ->setComment("Path ONNX Energy model");
    iDesc.add<double>("eid_min_cluster_energy", 1.0);
    iDesc.add<int>("eid_n_layers", 50);
    iDesc.add<int>("eid_n_clusters", 10);
  }
}
// Define this as a plug-in
//DEFINE_EDM_PLUGIN(TracksterInferenceAlgoFactory, ticl::TracksterInferenceByDNN, "TracksterInferenceByDNN");
