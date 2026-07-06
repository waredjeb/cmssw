// Author: Felice Pantaleo, felice.pantaleo@cern.ch 08/2024

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"
#include "DataFormats/Provenance/interface/ProductID.h"

class AllLayerClusterToTracksterAssociatorsProducer : public edm::global::EDProducer<> {
public:
  explicit AllLayerClusterToTracksterAssociatorsProducer(const edm::ParameterSet&);
  ~AllLayerClusterToTracksterAssociatorsProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  edm::EDGetTokenT<std::vector<reco::CaloCluster>> layerClustersToken_;
  // Optional per-trackster-collection layer clusters (parallel to
  // tracksterCollectionTokens_). Empty => use the single layer_clusters.
  std::vector<edm::EDGetTokenT<std::vector<reco::CaloCluster>>> layerClustersByCollectionTokens_;
  std::vector<std::pair<std::string, edm::EDGetTokenT<std::vector<ticl::Trackster>>>> tracksterCollectionTokens_;
};

AllLayerClusterToTracksterAssociatorsProducer::AllLayerClusterToTracksterAssociatorsProducer(
    const edm::ParameterSet& pset)
    : layerClustersToken_(
          consumes<std::vector<reco::CaloCluster>>(pset.getParameter<edm::InputTag>("layer_clusters"))) {
  const auto& tracksterCollections = pset.getParameter<std::vector<edm::InputTag>>("tracksterCollections");
  for (const auto& tag : tracksterCollections) {
    std::string label = tag.label();
    if (!tag.instance().empty()) {
      label += tag.instance();
    }
    tracksterCollectionTokens_.emplace_back(label, consumes<std::vector<ticl::Trackster>>(tag));
  }

  // Optional per-collection layer clusters, one-to-one with tracksterCollections.
  const auto& layerClustersByCollection =
      pset.getParameter<std::vector<edm::InputTag>>("layerClustersByCollection");
  if (!layerClustersByCollection.empty()) {
    if (layerClustersByCollection.size() != tracksterCollections.size()) {
      throw cms::Exception("Configuration")
          << "AllLayerClusterToTracksterAssociatorsProducer: 'layerClustersByCollection' has "
          << layerClustersByCollection.size() << " entries but 'tracksterCollections' has "
          << tracksterCollections.size() << ". They must match one-to-one.";
    }
    for (const auto& tag : layerClustersByCollection) {
      layerClustersByCollectionTokens_.emplace_back(consumes<std::vector<reco::CaloCluster>>(tag));
    }
  }

  // Produce separate association maps for each trackster collection using the trackster label
  for (const auto& tracksterToken : tracksterCollectionTokens_) {
    produces<
        ticl::AssociationMap<ticl::mapWithSharedEnergy, std::vector<reco::CaloCluster>, std::vector<ticl::Trackster>>>(
        tracksterToken.first);
  }
}

void AllLayerClusterToTracksterAssociatorsProducer::produce(edm::StreamID,
                                                            edm::Event& iEvent,
                                                            const edm::EventSetup&) const {
  using namespace edm;

  // Retrieve layer clusters with protection
  const auto& layerClustersHandle = iEvent.getHandle(layerClustersToken_);

  // If layer clusters are missing, produce empty maps and return
  if (!layerClustersHandle.isValid()) {
    edm::LogWarning("MissingInput") << "Layer clusters collection not found. Producing empty maps.";
    for (const auto& tracksterToken : tracksterCollectionTokens_) {
      iEvent.put(std::make_unique<ticl::AssociationMap<ticl::mapWithSharedEnergy,
                                                       std::vector<reco::CaloCluster>,
                                                       std::vector<ticl::Trackster>>>(),
                 tracksterToken.first);
    }
    return;
  }

  for (unsigned int ic = 0; ic < tracksterCollectionTokens_.size(); ++ic) {
    const auto& tracksterToken = tracksterCollectionTokens_[ic];
    const auto& trackstersHandle = iEvent.getHandle(tracksterToken.second);
    // If tracksters collection is missing, produce empty map and continue
    if (!trackstersHandle.isValid()) {
      LogDebug("MissingInput") << "Tracksters collection '" << tracksterToken.first << "' not found.";
      iEvent.put(std::make_unique<ticl::AssociationMap<ticl::mapWithSharedEnergy,
                                                       std::vector<reco::CaloCluster>,
                                                       std::vector<ticl::Trackster>>>(),
                 tracksterToken.first);
      continue;
    }

    // Resolve this collection's vertices against its own layer clusters when a
    // per-collection collection was configured, otherwise the shared one. The
    // association map is keyed by this layer-cluster collection.
    auto lcHandle = layerClustersHandle;
    if (!layerClustersByCollectionTokens_.empty()) {
      lcHandle = iEvent.getHandle(layerClustersByCollectionTokens_[ic]);
      if (!lcHandle.isValid()) {
        edm::LogWarning("MissingInput")
            << "Per-collection layer clusters for '" << tracksterToken.first << "' not found.";
        iEvent.put(std::make_unique<ticl::AssociationMap<ticl::mapWithSharedEnergy,
                                                         std::vector<reco::CaloCluster>,
                                                         std::vector<ticl::Trackster>>>(),
                   tracksterToken.first);
        continue;
      }
    }

    // Create association map
    auto lcToTracksterMap = std::make_unique<
        ticl::AssociationMap<ticl::mapWithSharedEnergy, std::vector<reco::CaloCluster>, std::vector<ticl::Trackster>>>(
        lcHandle, trackstersHandle, iEvent);

    // Loop over tracksters
    for (unsigned int tracksterId = 0; tracksterId < trackstersHandle->size(); ++tracksterId) {
      const auto& trackster = (*trackstersHandle)[tracksterId];
      // Loop over vertices in trackster
      for (unsigned int i = 0; i < trackster.vertices().size(); ++i) {
        // Get layerCluster
        const auto& lc = (*lcHandle)[trackster.vertices()[i]];
        float sharedEnergy = lc.energy() / trackster.vertex_multiplicity()[i];
        edm::Ref<std::vector<reco::CaloCluster>> lcRef(lcHandle, trackster.vertices()[i]);
        edm::Ref<std::vector<ticl::Trackster>> tracksterRef(trackstersHandle, tracksterId);
        lcToTracksterMap->insert(lcRef, tracksterRef, sharedEnergy);
      }
    }

    iEvent.put(std::move(lcToTracksterMap), tracksterToken.first);
  }
}

void AllLayerClusterToTracksterAssociatorsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::vector<edm::InputTag>>("tracksterCollections",
                                       {edm::InputTag("ticlTrackstersCLUE3DHigh"),
                                        edm::InputTag("ticlTrackstersLinks"),
                                        edm::InputTag("ticlCandidate")});
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalMergeLayerClusters"));
  // Optional per-trackster-collection layer clusters (parallel to
  // tracksterCollections). Empty by default => use the single 'layer_clusters'.
  desc.add<std::vector<edm::InputTag>>("layerClustersByCollection", {});
  descriptions.add("AllLayerClusterToTracksterAssociatorsProducer", desc);
}

// Define this as a plug-in
DEFINE_FWK_MODULE(AllLayerClusterToTracksterAssociatorsProducer);
