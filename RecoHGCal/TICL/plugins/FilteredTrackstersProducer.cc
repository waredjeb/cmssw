// Author: Felice Pantaleo, Marco Rovere - felice.pantaleo@cern.ch, marco.rovere@cern.ch
// Date: 09/2018

// user include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "TracksterFilterFactory.h"
#include "TracksterFilterBase.h"

#include <string>

class FilteredTrackstersProducer : public edm::stream::EDProducer<> {
public:
  FilteredTrackstersProducer(const edm::ParameterSet&);
  ~FilteredTrackstersProducer() override {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_token_;
  edm::EDGetTokenT<std::vector<float>> trackstersMask_token_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometry_token_;
  std::string tracksterFilter_;
  std::string iteration_label_;
  std::unique_ptr<const ticl::TracksterFilterBase> theFilter_;
  hgcal::RecHitTools rhtools_;
};

DEFINE_FWK_MODULE(FilteredTrackstersProducer);

FilteredTrackstersProducer::FilteredTrackstersProducer(const edm::ParameterSet& ps) {
  clusters_token_ = consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("LayerClusters"));
  tracksters_token_ = consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("Tracksters"));
  trackstersMask_token_ = consumes<std::vector<float>>(ps.getParameter<edm::InputTag>("TrackstersInputMask"));
  caloGeometry_token_ = esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>();
  tracksterFilter_ = ps.getParameter<std::string>("tracksterFilter");
  theFilter_ = TracksterFilterFactory::get()->create(tracksterFilter_, ps);
  iteration_label_ = ps.getParameter<std::string>("iteration_label");
  produces<std::vector<float>>(iteration_label_);
}

void FilteredTrackstersProducer::beginRun(edm::Run const&, edm::EventSetup const& es) {
  edm::ESHandle<CaloGeometry> geom = es.getHandle(caloGeometry_token_);
  rhtools_.setGeometry(*geom);
}

void FilteredTrackstersProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("LayerClusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("Tracksters", edm::InputTag("ticlTrackstersCLUE3DHigh"));
  desc.add<edm::InputTag>("TrackstersInputMask",
                          edm::InputTag("ticlTrackstersCLUE3DHigh", "tracksterMask"));
  desc.add<std::string>("iteration_label", "iterationLabelGoesHere");
  desc.add<std::string>("tracksterFilter", "TracksterFilterByPDGID");
  desc.add<std::vector<int>>(
      "algo_number",
      {reco::CaloCluster::hgcal_em, reco::CaloCluster::hgcal_had, reco::CaloCluster::hgcal_scintillator});  // 6,7,8
  desc.add<int>("min_cluster_size", 0);
  desc.add<int>("max_cluster_size", 9999);
  desc.add<bool>("filterEM", true);
  desc.add<double>("threshold", true);
  desc.add<int>("min_layerId", 0);
  desc.add<int>("max_layerId", 9999);

  descriptions.add("filteredTrackstersProducer", desc);
}

void FilteredTrackstersProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
  edm::Handle<std::vector<reco::CaloCluster>> clusterHandle;
  edm::Handle<std::vector<ticl::Trackster>> tracksterHandle;
  edm::Handle<std::vector<float>> inputTrackstersMaskHandle;
  evt.getByToken(clusters_token_, clusterHandle);
  evt.getByToken(tracksters_token_, tracksterHandle);
  evt.getByToken(trackstersMask_token_, inputTrackstersMaskHandle);

  // Protection against missing input collections
  if (!tracksterHandle.isValid() || !inputTrackstersMaskHandle.isValid()) {
    edm::LogWarning("FilteredTrackstersProducer") << "Missing input collections. Producing an empty mask.";

    // Produce an empty mask and exit
    auto emptyMask = std::make_unique<std::vector<float>>();
    evt.put(std::move(emptyMask), iteration_label_);
    return;
  }

  const auto& inputTrackstersMask = *inputTrackstersMaskHandle;

  // Transfer input mask in output
  auto trackstersMask = std::make_unique<std::vector<float>>(inputTrackstersMask);

  const auto& tracksters = *tracksterHandle;
  const auto& layerClusters = *clusterHandle;
  if (theFilter_) {
    theFilter_->filter(tracksters, layerClusters, *trackstersMask, rhtools_);
  }

  evt.put(std::move(trackstersMask), iteration_label_);
}
