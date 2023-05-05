#ifndef ENERGYREGRESSIONANDIDMODEL_H_
#define ENERGYREGRESSIONANDIDMODEL_H_

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

using namespace ticl;

class EnergyRegressionAndIDModel {
public:
  EnergyRegressionAndIDModel(const std::string eidInputName,
                             const std::string eidOutputNameEnergy,
                             const std::string eidOutputNameId,
                             const float eidMinClusterEnergy,
                             const int eidNLayers,
                             const int eidNClusters,
                             const tensorflow::Session *eidSession,
                             hgcal::RecHitTools &rhtools)
      : eidInputName_(eidInputName),
        eidOutputNameEnergy_(eidOutputNameEnergy),
        eidOutputNameId_(eidOutputNameId),
        eidMinClusterEnergy_(eidMinClusterEnergy),
        eidNLayers_(eidNLayers),
        eidNClusters_(eidNClusters),
        eidSession_(eidSession),
        rhtools_(rhtools){};

  ~EnergyRegressionAndIDModel() = default;
  void energyRegressionAndID(const std::vector<reco::CaloCluster> &layerClusters,
                             std::vector<Trackster> &tracksters) const {
    // Energy regression and particle identification strategy:
    //
    // 1. Set default values for regressed energy and particle id for each trackster.
    // 2. Store indices of tracksters whose total sum of cluster energies is above the
    //    eidMinClusterEnergy_ (GeV) treshold. Inference is not applied for soft tracksters.
    // 3. When no trackster passes the selection, return.
    // 4. Create input and output tensors. The batch dimension is determined by the number of
    //    selected tracksters.
    // 5. Fill input tensors with layer cluster features. Per layer, clusters are ordered descending
    //    by energy. Given that tensor data is contiguous in memory, we can use pointer arithmetic to
    //    fill values, even with batching.
    // 6. Zero-fill features for empty clusters in each layer.
    // 7. Batched inference.
    // 8. Assign the regressed energy and id probabilities to each trackster.
    //
    // Indices used throughout this method:
    // i -> batch element / trackster
    // j -> layer
    // k -> cluster
    // l -> feature

    // do nothing when no trackster passes the selection (3)
    int batchSize = (int)tracksters.size();
    if (batchSize == 0) {
      return;
    }

    for (auto &t : tracksters) {
      t.setRegressedEnergy(0.f);
      t.zeroProbabilities();
    }

    // create input and output tensors (4)
    tensorflow::TensorShape shape({batchSize, eidNLayers_, eidNClusters_, eidNFeatures_});
    tensorflow::Tensor input(tensorflow::DT_FLOAT, shape);
    tensorflow::NamedTensorList inputList = {{eidInputName_, input}};
    static constexpr int inputDimension = 4;

    std::vector<tensorflow::Tensor> outputs;
    std::vector<std::string> outputNames;
    if (!eidOutputNameEnergy_.empty()) {
      outputNames.push_back(eidOutputNameEnergy_);
    }
    if (!eidOutputNameId_.empty()) {
      outputNames.push_back(eidOutputNameId_);
    }

    // fill input tensor (5)
    for (int i = 0; i < batchSize; i++) {
      const Trackster &trackster = tracksters[i];

      // per layer, we only consider the first eidNClusters_ clusters in terms of
      // energy, so in order to avoid creating large / nested structures to do
      // the sorting for an unknown number of total clusters, create a sorted
      // list of layer cluster indices to keep track of the filled clusters
      std::vector<int> clusterIndices(trackster.vertices().size());
      for (int k = 0; k < (int)trackster.vertices().size(); k++) {
        clusterIndices[k] = k;
      }
      sort(clusterIndices.begin(), clusterIndices.end(), [&layerClusters, &trackster](const int &a, const int &b) {
        return layerClusters[trackster.vertices(a)].energy() > layerClusters[trackster.vertices(b)].energy();
      });

      // keep track of the number of seen clusters per layer
      std::vector<int> seenClusters(eidNLayers_);

      // loop through clusters by descending energy
      for (const int &k : clusterIndices) {
        // get features per layer and cluster and store the values directly in the input tensor
        const reco::CaloCluster &cluster = layerClusters[trackster.vertices(k)];
        int j = rhtools_.getLayerWithOffset(cluster.hitsAndFractions()[0].first) - 1;
        if (j < eidNLayers_ && seenClusters[j] < eidNClusters_) {
          // get the pointer to the first feature value for the current batch, layer and cluster
          float *features = &input.tensor<float, inputDimension>()(i, j, seenClusters[j], 0);

          // fill features
          *(features++) = float(cluster.energy() / float(trackster.vertex_multiplicity(k)));
          *(features++) = float(std::abs(cluster.eta()));
          *(features) = float(cluster.phi());

          // increment seen clusters
          seenClusters[j]++;
        }
      }

      // zero-fill features of empty clusters in each layer (6)
      for (int j = 0; j < eidNLayers_; j++) {
        for (int k = seenClusters[j]; k < eidNClusters_; k++) {
          float *features = &input.tensor<float, inputDimension>()(i, j, k, 0);
          for (int l = 0; l < eidNFeatures_; l++) {
            *(features++) = 0.f;
          }
        }
      }
    }

    // run the inference (7)
    tensorflow::run(const_cast<tensorflow::Session *>(eidSession_), inputList, outputNames, &outputs);

    // store regressed energy per trackster (8)
    if (!eidOutputNameEnergy_.empty()) {
      // get the pointer to the energy tensor, dimension is batch x 1
      float *energy = outputs[0].flat<float>().data();

      for (int i = 0; i < batchSize; ++i) {
        float regressedEnergy =
            tracksters[i].raw_energy() > eidMinClusterEnergy_ ? energy[i] : tracksters[i].raw_energy();
        tracksters[i].setRegressedEnergy(regressedEnergy);
      }
    }

    // store id probabilities per trackster (8)
    if (!eidOutputNameId_.empty()) {
      // get the pointer to the id probability tensor, dimension is batch x id_probabilities.size()
      int probsIdx = !eidOutputNameEnergy_.empty();
      float *probs = outputs[probsIdx].flat<float>().data();
      int probsNumber = tracksters[0].id_probabilities().size();
      for (int i = 0; i < batchSize; ++i) {
        tracksters[i].setProbabilities(&probs[i * probsNumber]);
      }
    }
  }

private:
  const std::string eidInputName_;
  const std::string eidOutputNameEnergy_;
  const std::string eidOutputNameId_;
  const float eidMinClusterEnergy_;
  const int eidNLayers_;
  const int eidNClusters_;
  const tensorflow::Session *eidSession_;
  static constexpr int eidNFeatures_ = 3;
  hgcal::RecHitTools rhtools_;
};
#endif  // ENERGYREGRESSIONANDIDMODEL_H_
