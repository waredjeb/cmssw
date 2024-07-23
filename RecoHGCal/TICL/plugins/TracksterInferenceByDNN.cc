#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h" 
#include "RecoHGCal/TICL/interface/TracksterInferenceByDNN.h" 
#include "RecoHGCal/TICL/interface/TracksterInferenceAlgoFactory.h" 
#include "FWCore/ParameterSet/interface/ParameterSet.h" 
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "RecoHGCal/TICL/interface/PatternRecognitionAlgoBase.h" 
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h" 
#include "TrackstersPCA.h" 

namespace ticl {
  using namespace cms::Ort; // Use ONNXRuntime namespace

  // Constructor for TracksterInferenceByDNN
  TracksterInferenceByDNN::TracksterInferenceByDNN(const edm::ParameterSet& conf)
    : TracksterInferenceAlgoBase(conf),
      id_modelPath(conf.getParameter<edm::FileInPath>("onnxPIDModelPath").fullPath()), // Path to the PID model
      en_modelPath(conf.getParameter<edm::FileInPath>("onnxEnergyModelPath").fullPath()), // Path to the Energy model
      eidMinClusterEnergy_(conf.getParameter<double>("eid_min_cluster_energy")), // Minimum cluster energy
      eidNLayers_(conf.getParameter<int>("eid_n_layers")), // Number of layers
      eidNClusters_(conf.getParameter<int>("eid_n_clusters")) // Number of clusters
  {
    // Initialize ONNX Runtime sessions for PID and Energy models
    static std::unique_ptr<cms::Ort::ONNXRuntime> onnxPIDRuntimeInstance = std::make_unique<cms::Ort::ONNXRuntime>(id_modelPath.c_str());
    onnxPIDSession_ = onnxPIDRuntimeInstance.get();

    static std::unique_ptr<cms::Ort::ONNXRuntime> onnxEnergyRuntimeInstance = std::make_unique<cms::Ort::ONNXRuntime>(en_modelPath.c_str());
    onnxEnergySession_ = onnxEnergyRuntimeInstance.get();
  }

  // Method to process input data and prepare it for inference
  void TracksterInferenceByDNN::inputData(const std::vector<reco::CaloCluster>& layerClusters, std::vector<Trackster>& tracksters) {
    tracksterIndices.clear(); // Clear previous indices
    for (int i = 0; i < static_cast<int>(tracksters.size()); i++) {
      float sumClusterEnergy = 0.;
      for (const unsigned int& vertex : tracksters[i].vertices()) {
        sumClusterEnergy += static_cast<float>(layerClusters[vertex].energy());
        if (sumClusterEnergy >= eidMinClusterEnergy_) {
          tracksters[i].setRegressedEnergy(0.f); // Set regressed energy to 0
          tracksters[i].zeroProbabilities(); // Zero out probabilities
          tracksterIndices.push_back(i); // Add index to the list
          break;
        }
      }
    }

    // Prepare input shapes and data for inference
    batchSize = static_cast<int>(tracksterIndices.size());
    if (batchSize == 0) return; // Exit if no tracksters

    std::vector<int64_t> inputShape = {batchSize, eidNLayers_, eidNClusters_, eidNFeatures_};
    input_shapes = {inputShape};

    input_Data.clear();
    input_Data.emplace_back(batchSize * eidNLayers_ * eidNClusters_ * eidNFeatures_, 0);

    for (int i = 0; i < batchSize; i++) {
      const Trackster& trackster = tracksters[tracksterIndices[i]];

      // Prepare indices and sort clusters based on energy
      std::vector<int> clusterIndices(trackster.vertices().size());
      for (int k = 0; k < static_cast<int>(trackster.vertices().size()); k++) {
        clusterIndices[k] = k;
      }

      std::sort(clusterIndices.begin(), clusterIndices.end(), [&layerClusters, &trackster](const int& a, const int& b) {
        return layerClusters[trackster.vertices(a)].energy() > layerClusters[trackster.vertices(b)].energy();
      });

      std::vector<int> seenClusters(eidNLayers_, 0);

      // Fill input data with cluster information
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

  // Method to run inference and update tracksters
  void TracksterInferenceByDNN::runInference(std::vector<Trackster>& tracksters) {
    if (batchSize == 0) return; // Exit if no batch

    // Define input and output names for inference
    std::vector<std::string> inputNames = {"input"};
    std::vector<std::string> output_en  = {"enreg_output"};
    std::vector<std::string> output_id  = {"pid_output"};

    // Run energy model inference
    std::vector<float> energyOutputTensor = onnxEnergySession_->run(inputNames, input_Data, input_shapes, output_en, batchSize)[0];

    if (!output_en.empty()) {
      for (int i = 0; i < batchSize; i++) {
        float energy = energyOutputTensor[i];
        tracksters[tracksterIndices[i]].setRegressedEnergy(energy); // Update energy
      }
    }

    // Run PID model inference
    std::vector<float> pidOutputTensor = onnxPIDSession_->run(inputNames, input_Data, input_shapes, output_id, batchSize)[0];
    
    if (!output_id.empty()) {
      float* probs = pidOutputTensor.data();
      for (int i = 0; i < batchSize; i++) {
        tracksters[tracksterIndices[i]].setProbabilities(probs); // Update probabilities
        probs += tracksters[tracksterIndices[i]].id_probabilities().size(); // Move to next set of probabilities
      }
    }
  }

  // Method to fill parameter set description for configuration
  void TracksterInferenceByDNN::fillPSetDescription(edm::ParameterSetDescription& iDesc) {
    iDesc.add<int>("algo_verbosity", 0); 
    iDesc.add<edm::FileInPath>("onnxPIDModelPath", edm::FileInPath("RecoHGCal/TICL/data/onnx_models/id_v0.onnx"))
      ->setComment("Path to ONNX PID model"); 
    iDesc.add<edm::FileInPath>("onnxEnergyModelPath", edm::FileInPath("RecoHGCal/TICL/data/onnx_models/energy_v0.onnx"))
      ->setComment("Path to ONNX Energy model"); 
    iDesc.add<double>("eid_min_cluster_energy", 1.0); 
    iDesc.add<int>("eid_n_layers", 50);
    iDesc.add<int>("eid_n_clusters", 10); 
  }
}
