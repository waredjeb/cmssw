// Author: Wahid Redjeb - wahid.redjeb@cern.ch
// Date: 01/2026

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "TracksterFilterFactory.h"
#include "TracksterFilterBase.h"

#include <string>
#include <memory>

class TracksterFilterProducer : public edm::stream::EDProducer<> {
public:
  TracksterFilterProducer(const edm::ParameterSet&);
  ~TracksterFilterProducer() override {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_token_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> layer_clusters_token_;
  edm::EDGetTokenT<std::vector<float>> tracksters_mask_token_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometry_token_;
  std::string filter_type_;
  std::unique_ptr<ticl::TracksterFilterBase> filter_;
  hgcal::RecHitTools rhtools_;
};

DEFINE_FWK_MODULE(TracksterFilterProducer);

TracksterFilterProducer::TracksterFilterProducer(const edm::ParameterSet& ps) {
  tracksters_token_ = consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("tracksters"));
  layer_clusters_token_ =
      consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layerClusters"));
  tracksters_mask_token_ = consumes<std::vector<float>>(ps.getParameter<edm::InputTag>("trackstersMask"));
  geometry_token_ = esConsumes<CaloGeometry, CaloGeometryRecord>();

  filter_type_ = ps.getParameter<std::string>("filterType");
  const edm::ParameterSet& filterPSet = ps.getParameter<edm::ParameterSet>("filterParams");
  filter_ = std::unique_ptr<ticl::TracksterFilterBase>(TracksterFilterFactory::get()->create(filter_type_, filterPSet));

  produces<std::vector<float>>();
}

void TracksterFilterProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracksters", edm::InputTag("ticlTrackstersCLUE3DHigh"))
      ->setComment("Input Trackster collection");
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalMergeLayerClusters"))
      ->setComment("Input LayerCluster collection");
  desc.add<edm::InputTag>("trackstersMask", edm::InputTag("ticlTrackstersCLUE3DHigh", "tracksterMask"))
      ->setComment("Input Trackster mask");
  desc.add<std::string>("filterType", "TracksterFilterByPDGID")->setComment("Name of the filter plugin to use");

  edm::ParameterSetDescription filterParamsDesc;
  filterParamsDesc.add<bool>("keepHadronic", true)->setComment("Keep hadronic (true) or EM (false) tracksters");
  desc.add<edm::ParameterSetDescription>("filterParams", filterParamsDesc)
      ->setComment("Parameters for the specific filter");

  descriptions.add("tracksterFilterProducer", desc);
}

void TracksterFilterProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
  edm::ESHandle<CaloGeometry> geom = es.getHandle(geometry_token_);
  rhtools_.setGeometry(*geom);

  edm::Handle<std::vector<ticl::Trackster>> trackstersHandle;
  edm::Handle<std::vector<reco::CaloCluster>> layerClustersHandle;
  edm::Handle<std::vector<float>> trackstersMaskHandle;

  evt.getByToken(tracksters_token_, trackstersHandle);
  evt.getByToken(layer_clusters_token_, layerClustersHandle);
  evt.getByToken(tracksters_mask_token_, trackstersMaskHandle);

  if (!trackstersHandle.isValid() || !trackstersMaskHandle.isValid() || !layerClustersHandle.isValid()) {
    edm::LogWarning("TracksterFilterProducer") << "Missing input collections. Producing an empty mask.";
    auto emptyMask = std::make_unique<std::vector<float>>();
    evt.put(std::move(emptyMask));
    return;
  }

  const auto& tracksters = *trackstersHandle;
  const auto& layerClusters = *layerClustersHandle;
  const auto& inputMask = *trackstersMaskHandle;

  auto filteredMask = std::make_unique<std::vector<float>>(inputMask);

  // apply filter 
  if (filter_) {
    filter_->filter(tracksters, layerClusters, *filteredMask, rhtools_);
  }
  std::cout << "TracksterFilterProducer input collection size " << inputMask.size() << " Output size " << filteredMask->size() << std::endl; 
  evt.put(std::move(filteredMask));
}
