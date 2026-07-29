// Authors: Olivie Franklova - olivie.abigail.franklova@cern.ch
// Date: 03/2023
// @file merge layer clusters which were produce by HGCalLayerClusterProducer

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterHostCollection.h"
#include "DataFormats/TICL/interface/HitsAndFractionsHost.h"

#include "HeterogeneousCore/AlpakaInterface/interface/host.h"

#include <cassert>
#include <vector>

class MergeClusterProducer : public edm::stream::EDProducer<> {
public:
  /**
   * @brief Constructor with parameter settings - which can be changed in  ...todo.
   * Constructor will set all variables by input param ps.
   *
   * @param[in] ps parametr set to set variables
  */
  MergeClusterProducer(const edm::ParameterSet &);
  ~MergeClusterProducer() override {}
  /**
   * @brief Method fill description which will be used in pyhton file.
   *
   * @param[out] description to be fill
  */
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  /**
   * @brief Method will merge the producers and put them back to event
   *
   * @param[in, out] evt from get info and put result
   * @param[in] es to get event setup info
  */
  void produce(edm::Event &, const edm::EventSetup &) override;

private:
  // Each layer-cluster producer emits a cluster SoA and, under the same label,
  // the matching hits-and-fractions association map. The two are consumed from
  // the same InputTag: EDM resolves products by type.
  std::vector<edm::EDGetTokenT<reco::CaloClusterHostCollection>> tokens_;
  std::vector<edm::EDGetTokenT<ticl::HitsAndFractionsHost>> hits_tokens_;

  /**
   * @brief Append one input cluster SoA to the merged one
   *
   * @param[out] merged the collection the input is appended to
   * @param[in] input the collection to append
   * @param[in] start index in merged at which the input is written
  */
  static void mergeClusters(reco::CaloClusterHostCollection::View &merged,
                            const reco::CaloClusterHostCollection::ConstView &input,
                            int start);

  /**
   * @brief Append one input association map to the merged one
   *
   * The map is keyed by cluster index, so the keys are shifted by the same
   * amount as the clusters (clusterStart) and the offsets by the number of hits
   * already written (hitStart). Keeping these two in step with mergeClusters is
   * what makes merged[i] the hit list of merged cluster i.
  */
  static void mergeHitsAndFractions(ticl::HitsAndFractionsHost::View &merged,
                                    const ticl::HitsAndFractionsHost::ConstView &input,
                                    int clusterStart,
                                    int hitStart);
};

MergeClusterProducer::MergeClusterProducer(const edm::ParameterSet &ps) {
  std::vector<edm::InputTag> tags = ps.getParameter<std::vector<edm::InputTag>>("layerClusters");
  for (auto &tag : tags) {
    tokens_.push_back(consumes<reco::CaloClusterHostCollection>(tag));
    hits_tokens_.push_back(consumes<ticl::HitsAndFractionsHost>(tag));
  }

  produces<std::vector<float>>("InitialLayerClustersMask");
  produces<reco::CaloClusterHostCollection>();
  produces<ticl::HitsAndFractionsHost>();
}

void MergeClusterProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // hgcalMergeLayerClusters
  edm::ParameterSetDescription desc;
  //layer clusters
  desc.add<std::vector<edm::InputTag>>("layerClusters",
                                       {edm::InputTag("hgcalLayerClustersEE"),
                                        edm::InputTag("hgcalLayerClustersHSi"),
                                        edm::InputTag("hgcalLayerClustersHSci")});
  descriptions.add("hgcalMergeLayerClusters", desc);
}

void MergeClusterProducer::produce(edm::Event &evt, const edm::EventSetup &es) {
  std::vector<edm::Handle<reco::CaloClusterHostCollection>> clusterHandles;
  std::vector<edm::Handle<ticl::HitsAndFractionsHost>> hitsHandles;
  int totalClusters = 0;
  int totalHits = 0;

  for (size_t i = 0; i < tokens_.size(); ++i) {
    auto const &clusterHandle = evt.getHandle(tokens_[i]);
    auto const &hitsHandle = evt.getHandle(hits_tokens_[i]);
    totalClusters += clusterHandle->view().position().metadata().size();
    totalHits += hitsHandle->const_view().size();
    clusterHandles.push_back(clusterHandle);
    hitsHandles.push_back(hitsHandle);
  }

  auto merged = std::make_unique<reco::CaloClusterHostCollection>(
      cms::alpakatools::host(), totalClusters, totalClusters, totalClusters, totalClusters);
  auto mergedHits = std::make_unique<ticl::HitsAndFractionsHost>(cms::alpakatools::host(), totalHits, totalClusters);
  auto merged_v = merged->view();
  auto mergedHits_v = mergedHits->view();

  int clusterStart = 0;
  int hitStart = 0;
  for (size_t i = 0; i < clusterHandles.size(); ++i) {
    auto const input_v = clusterHandles[i]->const_view();
    auto const inputHits_v = hitsHandles[i]->const_view();
    mergeClusters(merged_v, input_v, clusterStart);
    mergeHitsAndFractions(mergedHits_v, inputHits_v, clusterStart, hitStart);
    clusterStart += input_v.position().metadata().size();
    hitStart += inputHits_v.size();
  }
  // CSR terminator: count(lastCluster) needs offsets[nClusters].
  mergedHits_v.offsets()[totalClusters].keys_offsets() = hitStart;
  // The map is keyed by position in the merged cluster collection, so a
  // mismatch here would silently mis-associate hits rather than crash.
  assert(clusterStart == totalClusters);
  assert(hitStart == totalHits);

  //create layer cluster mask
  auto layerClustersMask = std::make_unique<std::vector<float>>(totalClusters, 1.0);

  evt.put(std::move(merged));
  evt.put(std::move(mergedHits));
  evt.put(std::move(layerClustersMask), "InitialLayerClustersMask");
}

void MergeClusterProducer::mergeClusters(reco::CaloClusterHostCollection::View &merged,
                                         const reco::CaloClusterHostCollection::ConstView &input,
                                         int start) {
  for (int idx = 0; idx < input.position().metadata().size(); ++idx) {
    const auto cumulative_index = idx + start;
    merged.position()[cumulative_index].x() = input.position()[idx].x();
    merged.position()[cumulative_index].y() = input.position()[idx].y();
    merged.position()[cumulative_index].z() = input.position()[idx].z();
    merged.position()[cumulative_index].layer() = input.position()[idx].layer();
    merged.position()[cumulative_index].cells() = input.position()[idx].cells();
    merged.energy()[cumulative_index].energy() = input.energy()[idx].energy();
    merged.energy()[cumulative_index].correctedEnergy() = input.energy()[idx].correctedEnergy();
    merged.energy()[cumulative_index].correctedEnergyUncertainty() = input.energy()[idx].correctedEnergyUncertainty();
    merged.indexes()[cumulative_index].caloID() = input.indexes()[idx].caloID();
    merged.indexes()[cumulative_index].algoID() = input.indexes()[idx].algoID();
    merged.indexes()[cumulative_index].seedID() = input.indexes()[idx].seedID();
    merged.indexes()[cumulative_index].flags() = input.indexes()[idx].flags();
    merged.timing()[cumulative_index].time() = input.timing()[idx].time();
    merged.timing()[cumulative_index].timeError() = input.timing()[idx].timeError();
  }
}

void MergeClusterProducer::mergeHitsAndFractions(ticl::HitsAndFractionsHost::View &merged,
                                                 const ticl::HitsAndFractionsHost::ConstView &input,
                                                 int clusterStart,
                                                 int hitStart) {
  for (int key = 0; key < input.keys(); ++key) {
    merged.offsets()[key + clusterStart].keys_offsets() = input.offsets()[key].keys_offsets() + hitStart;
  }
  for (int idx = 0; idx < input.size(); ++idx) {
    merged.content()[idx + hitStart].values() = input.content()[idx].values();
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MergeClusterProducer);
