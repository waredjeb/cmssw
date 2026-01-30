// Author: Felice Pantaleo, Wahid Redjeb (CERN) - felice.pantaleo@cern.ch, wahid.redjeb@cern.ch
// Date: 12/2023
#include <memory>  // unique_ptr
#include "DataFormats/Common/interface/MultiSpan.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/Common/interface/OrphanHandle.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"

#include "PhysicsTools/TensorFlow/interface/TfGraphRecord.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "PhysicsTools/TensorFlow/interface/TfGraphDefWrapper.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

#include "RecoHGCal/TICL/interface/TracksterLinkingAlgoBase.h"
#include "RecoHGCal/TICL/plugins/TracksterLinkingPluginFactory.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "RecoHGCal/TICL/interface/TracksterInferenceAlgoFactory.h"

#include "TrackstersPCA.h"

using namespace ticl;
using cms::Ort::ONNXRuntime;

class TracksterLinksProducer : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>> {
public:
  explicit TracksterLinksProducer(const edm::ParameterSet &ps, const ONNXRuntime *);
  ~TracksterLinksProducer() override {};
  void produce(edm::Event &, const edm::EventSetup &) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  void beginRun(edm::Run const &iEvent, edm::EventSetup const &es) override;
  static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet &iConfig);
  static void globalEndJob(const ONNXRuntime *);

private:
  void printTrackstersDebug(const std::vector<Trackster> &, const char *label) const;
  void dumpTrackster(const Trackster &) const;

  std::unique_ptr<TracksterLinkingAlgoBase> linkingAlgo_;
  std::string algoType_;

  std::vector<edm::EDGetTokenT<std::vector<Trackster>>> tracksters_tokens_;
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;
  const edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> clustersTime_token_;

  const bool regressionAndPid_;
  std::unique_ptr<TracksterInferenceAlgoBase> inferenceAlgo_;

  std::vector<edm::EDGetTokenT<std::vector<float>>> original_masks_tokens_;
  std::vector<edm::EDGetTokenT<std::vector<float>>> original_tracksters_masks_tokens_;
  std::vector<edm::EDGetTokenT<std::vector<float>>> tracksters_masks_tokens_;
  std::vector<std::string> tracksters_collections_labels_;

  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometry_token_;
  const std::string detector_;
  const std::string propName_;

  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfield_token_;
  const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagator_token_;
  const HGCalDDDConstants *hgcons_;
  hgcal::RecHitTools rhtools_;
  edm::ESGetToken<HGCalDDDConstants, IdealGeometryRecord> hdc_token_;
};

TracksterLinksProducer::TracksterLinksProducer(const edm::ParameterSet &ps, const ONNXRuntime *onnxRuntime)
    : algoType_(ps.getParameter<edm::ParameterSet>("linkingPSet").getParameter<std::string>("type")),
      clusters_token_(consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_clusters"))),
      clustersTime_token_(
          consumes<edm::ValueMap<std::pair<float, float>>>(ps.getParameter<edm::InputTag>("layer_clustersTime"))),
      regressionAndPid_(ps.getParameter<bool>("regressionAndPid")),
      geometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>()),
      detector_(ps.getParameter<std::string>("detector")),
      propName_(ps.getParameter<std::string>("propagator")),
      bfield_token_(esConsumes<MagneticField, IdealMagneticFieldRecord, edm::Transition::BeginRun>()),
      propagator_token_(
          esConsumes<Propagator, TrackingComponentsRecord, edm::Transition::BeginRun>(edm::ESInputTag("", propName_))) {
  // Loop over the edm::VInputTag and append the token to tracksters_tokens_
  for (auto const &tag : ps.getParameter<std::vector<edm::InputTag>>("tracksters_collections")) {
    tracksters_tokens_.emplace_back(consumes<std::vector<Trackster>>(tag));
    tracksters_collections_labels_.emplace_back(tag.label());
  }
  //Loop over the edm::VInputTag of masks and append the token to original_masks_tokens_
  for (auto const &tag : ps.getParameter<std::vector<edm::InputTag>>("original_masks")) {
    original_masks_tokens_.emplace_back(consumes<std::vector<float>>(tag));
  }
  //Loop over the edm::VInputTag of original trackster masks and append the token to original_tracksters_masks_tokens_
  for (auto const &tag : ps.getParameter<std::vector<edm::InputTag>>("original_trackstersMasks")) {
    original_tracksters_masks_tokens_.emplace_back(consumes<std::vector<float>>(tag));
  }
  //Loop over the edm::VInputTag of filtered trackster masks and append the token to tracksters_masks_tokens_
  for (auto const &tag : ps.getParameter<std::vector<edm::InputTag>>("trackstersMasks")) {
    tracksters_masks_tokens_.emplace_back(consumes<std::vector<float>>(tag));
  }

  // Validation: ensure trackster masks match trackster collections
  auto const &tagTracksters = ps.getParameter<std::vector<edm::InputTag>>("tracksters_collections");
  auto const &tagMasks = ps.getParameter<std::vector<edm::InputTag>>("trackstersMasks");
  auto const &tagOriginalMasks = ps.getParameter<std::vector<edm::InputTag>>("original_trackstersMasks");
  assert(tagTracksters.size() == tagMasks.size());
  assert(tagTracksters.size() == tagOriginalMasks.size());
  if (tracksters_masks_tokens_.size() != tracksters_tokens_.size()) {
    throw cms::Exception("Configuration")
        << "Number of trackstersMasks (" << tracksters_masks_tokens_.size()
        << ") must match number of tracksters_collections (" << tracksters_tokens_.size() << ")";
  }
  if (original_tracksters_masks_tokens_.size() != tracksters_tokens_.size()) {
    throw cms::Exception("Configuration")
        << "Number of original_trackstersMasks (" << original_tracksters_masks_tokens_.size()
        << ") must match number of tracksters_collections (" << tracksters_tokens_.size() << ")";
  }
  // Initialize inference algorithm using the factory
  std::string inferencePlugin = ps.getParameter<std::string>("inferenceAlgo");
  edm::ParameterSet inferencePSet = ps.getParameter<edm::ParameterSet>("pluginInferenceAlgo" + inferencePlugin);
  inferenceAlgo_ = std::unique_ptr<TracksterInferenceAlgoBase>(
      TracksterInferenceAlgoFactory::get()->create(inferencePlugin, inferencePSet));

  // new trackster collection after linking
  produces<std::vector<Trackster>>();

  // links
  produces<std::vector<std::vector<unsigned int>>>();
  produces<std::vector<std::vector<unsigned int>>>("linkedTracksterIdToInputTracksterId");
  // layerClusters Mask
  produces<std::vector<float>>();
  // trackster mask for output linked tracksters
  produces<std::vector<float>>("tracksterMask");
  // updated trackster masks for each input collection
  for (const auto &label : tracksters_collections_labels_) {
    produces<std::vector<float>>("tracksterMask" + label);
  }

  auto linkingPSet = ps.getParameter<edm::ParameterSet>("linkingPSet");

  if (algoType_ == "Skeletons") {
    std::string detectorName_ = (detector_ == "HFNose") ? "HGCalHFNoseSensitive" : "HGCalEESensitive";
    hdc_token_ = esConsumes<HGCalDDDConstants, IdealGeometryRecord, edm::Transition::BeginRun>(
        edm::ESInputTag("", detectorName_));
  }

  linkingAlgo_ = TracksterLinkingPluginFactory::get()->create(algoType_, linkingPSet, consumesCollector(), onnxRuntime);
}

std::unique_ptr<ONNXRuntime> TracksterLinksProducer::initializeGlobalCache(const edm::ParameterSet &iConfig) {
  auto const &pluginPset = iConfig.getParameter<edm::ParameterSet>("linkingPSet");
  if (pluginPset.exists("onnxModelPath"))
    return std::make_unique<ONNXRuntime>(pluginPset.getParameter<edm::FileInPath>("onnxModelPath").fullPath());
  else
    return std::unique_ptr<ONNXRuntime>(nullptr);
}

void TracksterLinksProducer::globalEndJob(const ONNXRuntime *) {}

void TracksterLinksProducer::beginRun(edm::Run const &iEvent, edm::EventSetup const &es) {
  if (algoType_ == "Skeletons") {
    edm::ESHandle<HGCalDDDConstants> hdc = es.getHandle(hdc_token_);
    hgcons_ = hdc.product();
  }

  edm::ESHandle<CaloGeometry> geom = es.getHandle(geometry_token_);
  rhtools_.setGeometry(*geom);

  edm::ESHandle<MagneticField> bfield = es.getHandle(bfield_token_);
  edm::ESHandle<Propagator> propagator = es.getHandle(propagator_token_);

  linkingAlgo_->initialize(hgcons_, rhtools_, bfield, propagator);
};

void TracksterLinksProducer::dumpTrackster(const Trackster &t) const {
  auto e_over_h = (t.raw_em_pt() / ((t.raw_pt() - t.raw_em_pt()) != 0. ? (t.raw_pt() - t.raw_em_pt()) : 1.));
  LogDebug("TracksterLinksProducer")
      << "\nTrackster raw_pt: " << t.raw_pt() << " raw_em_pt: " << t.raw_em_pt() << " eoh: " << e_over_h
      << " barycenter: " << t.barycenter() << " eta,phi (baricenter): " << t.barycenter().eta() << ", "
      << t.barycenter().phi() << " eta,phi (eigen): " << t.eigenvectors(0).eta() << ", " << t.eigenvectors(0).phi()
      << " pt(eigen): " << std::sqrt(t.eigenvectors(0).Unit().perp2()) * t.raw_energy() << " seedID: " << t.seedID()
      << " seedIndex: " << t.seedIndex() << " size: " << t.vertices().size() << " average usage: "
      << (std::accumulate(std::begin(t.vertex_multiplicity()), std::end(t.vertex_multiplicity()), 0.) /
          (float)t.vertex_multiplicity().size())
      << " raw_energy: " << t.raw_energy() << " regressed energy: " << t.regressed_energy()
      << " probs(ga/e/mu/np/cp/nh/am/unk): ";
  for (auto const &p : t.id_probabilities()) {
    LogDebug("TracksterLinksProducer") << std::fixed << p << " ";
  }
  LogDebug("TracksterLinksProducer") << " sigmas: ";
  for (auto const &s : t.sigmas()) {
    LogDebug("TracksterLinksProducer") << s << " ";
  }
  LogDebug("TracksterLinksProducer") << std::endl;
}

void TracksterLinksProducer::produce(edm::Event &evt, const edm::EventSetup &es) {
  linkingAlgo_->setEvent(evt, es);
  auto resultTracksters = std::make_unique<std::vector<Trackster>>();

  auto linkedResultTracksters = std::make_unique<std::vector<std::vector<unsigned int>>>();

  const auto &layerClusters = evt.get(clusters_token_);
  const auto &layerClustersTimes = evt.get(clustersTime_token_);

  // loop over the original_masks_tokens_ and get the original masks collections and multiply them
  // to get the global mask
  std::vector<float> original_global_mask(layerClusters.size(), 1.f);
  for (unsigned int i = 0; i < original_masks_tokens_.size(); ++i) {
    const auto &tmp_mask = evt.get(original_masks_tokens_[i]);
    for (unsigned int j = 0; j < tmp_mask.size(); ++j) {
      original_global_mask[j] *= tmp_mask[j];
    }
  }

  auto resultMask = std::make_unique<std::vector<float>>(original_global_mask);

  std::vector<edm::Handle<std::vector<Trackster>>> tracksters_h(tracksters_tokens_.size());
  edm::MultiSpan<Trackster> trackstersManager;
  for (unsigned int i = 0; i < tracksters_tokens_.size(); ++i) {
    evt.getByToken(tracksters_tokens_[i], tracksters_h[i]);
    //Fill MultiSpan
    trackstersManager.add(*tracksters_h[i]);
    std::cout << "Adding trackster collection with size " << tracksters_h[i]->size() << std::endl;
  }
  std::cout << "Final Tracksters Manager size " << trackstersManager.size() << std::endl;

  // Get original trackster masks (all 1s)
  std::vector<std::vector<float>> originalTrackstersMasks;
  originalTrackstersMasks.reserve(original_tracksters_masks_tokens_.size());
  for (unsigned int i = 0; i < original_tracksters_masks_tokens_.size(); ++i) {
    const auto &originalMask = evt.get(original_tracksters_masks_tokens_[i]);
    originalTrackstersMasks.emplace_back(originalMask.begin(), originalMask.end());
    std::cout << "Adding Original Tracksters Mask collection with size " << originalMask.size() << std::endl;
  }

  // Get filtered input trackster masks and copy them for modification by linkTracksters
  std::vector<std::vector<float>> trackstersMasks;
  trackstersMasks.reserve(tracksters_masks_tokens_.size());
  size_t totalSizeMask = 0uz;
  for (unsigned int i = 0; i < tracksters_masks_tokens_.size(); ++i) {
    const auto &inputMask = evt.get(tracksters_masks_tokens_[i]);
    trackstersMasks.emplace_back(inputMask.begin(), inputMask.end());
    std::cout << "Adding Filtered Tracksters Mask collection with size " << inputMask.size() << std::endl;
    totalSizeMask += inputMask.size();
  }

  std::cout << "Final Tracksters mask size " << totalSizeMask << std::endl;

  for (size_t i{}; i < trackstersMasks.size(); i++) {
    assert(tracksters_h[i]->size() == trackstersMasks[i].size());
  }

  std::cout << "Trackster Links Producer - Original Mask " << algoType_ << std::endl;

  std::cout << "[";
  for (auto const &x : originalTrackstersMasks) {
    for (auto const i : x) {
      std::cout << i << ", ";
    }
  }
  std::cout << "]\n";
  std::cout << "Trackster Links Producer  - Input Mask " << algoType_ << std::endl;
  std::cout << "[";
  int j = 0;
  std::cout << " TrackstersMask Size " << trackstersMasks.size() << std::endl;
  for (auto const x : trackstersMasks) {
    std::cout << j << ", "; 
    for (auto const i : x) {
      std::cout << i << ", ";
    }
    j++;
  }
  std::cout << "]\n";

  // Save filtered masks before linking to determine which tracksters were used
  std::vector<std::vector<float>> filteredMasksBeforeLinking = trackstersMasks;

  // Linking
  const typename TracksterLinkingAlgoBase::Inputs input(evt, es, layerClusters, layerClustersTimes, trackstersManager);
  auto linkedTracksterIdToInputTracksterId = std::make_unique<std::vector<std::vector<unsigned int>>>();

  // LinkTracksters will produce a vector of vector of indices of tracksters that:
  // 1) are linked together if more than one
  // 2) are isolated if only one
  // Result tracksters contains the final version of the trackster collection
  // linkedTrackstersToInputTrackstersMap contains the mapping between the linked tracksters and the input tracksters
  // trackstersMasks will be modified to mark used tracksters as 0
  linkingAlgo_->linkTracksters(
      input, *resultTracksters, *linkedResultTracksters, *linkedTracksterIdToInputTracksterId, trackstersMasks);

  // Now we need to remove the tracksters that are not linked
  // We need to emplace_back in the resultTracksters only the tracksters that are linked

  for (auto const &resultTrackster : *resultTracksters) {
    for (auto const &clusterIndex : resultTrackster.vertices()) {
      (*resultMask)[clusterIndex] = 0.f;
    }
  }

  assignPCAtoTracksters(*resultTracksters,
                        layerClusters,
                        layerClustersTimes,
                        rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z(),
                        rhtools_,
                        true);

  if (regressionAndPid_) {
    // Run inference algorithm
    inferenceAlgo_->inputData(layerClusters, *resultTracksters, rhtools_);
    inferenceAlgo_->runInference(
        *resultTracksters);  //option to use "Linking" instead of "CLU3D"/"energyAndPid" instead of "PID"
  }
  std::cout << "Producing Tracksters Number : " << resultTracksters->size() << std::endl;
  std::cout << "Using CLUE3D Tracksters : \n";
  for (size_t i{}; i < resultTracksters->size(); i++){
    std::cout << "\t LinkedTrackster " << i << ": ";
    for (auto const iT : (*linkedTracksterIdToInputTracksterId)[i]){
      std::cout << (*linkedTracksterIdToInputTracksterId)[i][iT] << " "; 
    }
    std::cout << "\n";
  }

  // create output mask for linked tracksters 
  auto outputTrackstersMask = std::make_unique<std::vector<float>>(resultTracksters->size(), 1.f);

  evt.put(std::move(linkedResultTracksters));
  evt.put(std::move(resultMask));
  evt.put(std::move(resultTracksters));
  evt.put(std::move(linkedTracksterIdToInputTracksterId), "linkedTracksterIdToInputTracksterId");
  evt.put(std::move(outputTrackstersMask), "tracksterMask");



  // put updated trackster masks for each input collection
  // start from original masks and mark tracksters that were used
  for (unsigned int i = 0; i < originalTrackstersMasks.size(); ++i) {
    auto updatedMask = std::make_unique<std::vector<float>>(originalTrackstersMasks[i]);
    // mark tracksters that were used 
    for (unsigned int j = 0; j < trackstersMasks[i].size(); ++j) {
      // if trackster was available before (filtered mask = 1) and is now used (filtered mask = 0 after linking)
      if (filteredMasksBeforeLinking[i][j] > 0.f && trackstersMasks[i][j] == 0.f) {
        (*updatedMask)[j] = 0.f;
      }
    }
    std::cout << "Trackster Links Producer  - Output Updated Mask " << algoType_ << std::endl;
    std::cout << "[";
    for (auto const &i : *updatedMask) {
      std::cout << i << ", ";
    }
    std::cout << "]\n";
    evt.put(std::move(updatedMask), "tracksterMask" + tracksters_collections_labels_[i]);
  }
}

void TracksterLinksProducer::printTrackstersDebug(const std::vector<Trackster> &tracksters, const char *label) const {
  int counter = 0;
  LogDebug("TracksterLinksProducer").log([&](auto &log) {
    for (auto const &t : tracksters) {
      log << counter++ << " TracksterLinksProducer (" << label << ") obj barycenter: " << t.barycenter()
          << " eta,phi (baricenter): " << t.barycenter().eta() << ", " << t.barycenter().phi()
          << " eta,phi (eigen): " << t.eigenvectors(0).eta() << ", " << t.eigenvectors(0).phi()
          << " pt(eigen): " << std::sqrt(t.eigenvectors(0).Unit().perp2()) * t.raw_energy() << " seedID: " << t.seedID()
          << " seedIndex: " << t.seedIndex() << " size: " << t.vertices().size() << " average usage: "
          << (std::accumulate(std::begin(t.vertex_multiplicity()), std::end(t.vertex_multiplicity()), 0.) /
              (float)t.vertex_multiplicity().size())
          << " raw_energy: " << t.raw_energy() << " regressed energy: " << t.regressed_energy()
          << " probs(ga/e/mu/np/cp/nh/am/unk): ";
      for (auto const &p : t.id_probabilities()) {
        log << std::fixed << p << " ";
      }
      log << " sigmas: ";
      for (auto const &s : t.sigmas()) {
        log << s << " ";
      }
      log << "\n";
    }
  });
}

void TracksterLinksProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  edm::ParameterSetDescription linkingDesc;
  linkingDesc.addNode(edm::PluginDescription<TracksterLinkingPluginFactory>("type", "Skeletons", true));
  // Inference Plugins
  edm::ParameterSetDescription inferenceDesc;
  inferenceDesc.addNode(edm::PluginDescription<TracksterInferenceAlgoFactory>("type", "TracksterInferenceByDNN", true));
  desc.add<edm::ParameterSetDescription>("pluginInferenceAlgoTracksterInferenceByDNN", inferenceDesc);

  edm::ParameterSetDescription inferenceDescPFN;
  inferenceDescPFN.addNode(
      edm::PluginDescription<TracksterInferenceAlgoFactory>("type", "TracksterInferenceByPFN", true));
  desc.add<edm::ParameterSetDescription>("pluginInferenceAlgoTracksterInferenceByPFN", inferenceDescPFN);

  edm::ParameterSetDescription inferenceDescCNNv4;
  inferenceDescCNNv4.addNode(
      edm::PluginDescription<TracksterInferenceAlgoFactory>("type", "TracksterInferenceByCNNv4", true));
  desc.add<edm::ParameterSetDescription>("pluginInferenceAlgoTracksterInferenceByCNNv4", inferenceDescCNNv4);

  desc.add<edm::ParameterSetDescription>("linkingPSet", linkingDesc);
  desc.add<std::vector<edm::InputTag>>("tracksters_collections", {edm::InputTag("ticlTrackstersCLUE3DHigh")});
  desc.add<std::vector<edm::InputTag>>("original_masks",
                                       {edm::InputTag("hgcalMergeLayerClusters", "InitialLayerClustersMask")});
  desc.add<std::vector<edm::InputTag>>("original_trackstersMasks",
                                       {edm::InputTag("ticlTrackstersCLUE3DHigh", "tracksterMask")})
      ->setComment("Original trackster masks (all 1s) corresponding to each trackster collection");
  desc.add<std::vector<edm::InputTag>>("trackstersMasks", {edm::InputTag("filteredTrackstersCLUE3DHigh")})
      ->setComment("Filtered trackster masks corresponding to each trackster collection");
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("layer_clustersTime", edm::InputTag("hgcalMergeLayerClusters", "timeLayerCluster"));
  desc.add<bool>("regressionAndPid", false);
  desc.add<std::string>("detector", "HGCAL");
  desc.add<std::string>("propagator", "PropagatorWithMaterial");
  desc.add<std::string>("inferenceAlgo", "TracksterInferenceByPFN");
  descriptions.add("tracksterLinksProducer", desc);
}

DEFINE_FWK_MODULE(TracksterLinksProducer);
