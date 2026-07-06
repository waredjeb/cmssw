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
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/Common/interface/RefProdVector.h"
#include "DataFormats/Common/interface/MultiSpan.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"

template <typename HIT>
class AllHitToTracksterAssociatorsProducerT : public edm::global::EDProducer<> {
public:
  using multiCollectionT = edm::RefProdVector<std::vector<HIT>>;

  explicit AllHitToTracksterAssociatorsProducerT(const edm::ParameterSet&);
  ~AllHitToTracksterAssociatorsProducerT() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  std::vector<std::pair<std::string, edm::EDGetTokenT<std::vector<ticl::Trackster>>>> tracksterCollectionTokens_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> layerClustersToken_;
  // Optional per-trackster-collection layer clusters, parallel to
  // tracksterCollectionTokens_. Trackster collections that index into a
  // different layer-cluster collection (e.g. the L1-seeded tracksters, which
  // index into hltMergeLayerClustersL1Seeded) must resolve their vertices
  // against their own layer clusters. When empty, every collection falls back
  // to the single layerClustersToken_ (previous behaviour).
  std::vector<edm::EDGetTokenT<std::vector<reco::CaloCluster>>> layerClustersByCollectionTokens_;
  edm::EDGetTokenT<std::unordered_map<DetId, const unsigned int>> hitMapToken_;
  edm::EDGetTokenT<multiCollectionT> hitsToken_;
};

template <typename HIT>
AllHitToTracksterAssociatorsProducerT<HIT>::AllHitToTracksterAssociatorsProducerT(const edm::ParameterSet& pset)
    : layerClustersToken_(consumes<std::vector<reco::CaloCluster>>(pset.getParameter<edm::InputTag>("layerClusters"))),
      hitMapToken_(
          consumes<std::unordered_map<DetId, const unsigned int>>(pset.getParameter<edm::InputTag>("hitMapTag"))),
      hitsToken_(consumes<multiCollectionT>(pset.getParameter<edm::InputTag>("hits"))) {
  const auto& tracksterCollections = pset.getParameter<std::vector<edm::InputTag>>("tracksterCollections");
  for (const auto& tag : tracksterCollections) {
    tracksterCollectionTokens_.emplace_back(tag.label() + tag.instance(), consumes<std::vector<ticl::Trackster>>(tag));
  }

  // Optional per-collection layer clusters. If provided it must have the same
  // size as tracksterCollections; otherwise it is left empty and the single
  // layerClusters collection is used for every trackster collection.
  const auto& layerClustersByCollection =
      pset.getParameter<std::vector<edm::InputTag>>("layerClustersByCollection");
  if (!layerClustersByCollection.empty()) {
    if (layerClustersByCollection.size() != tracksterCollections.size()) {
      throw cms::Exception("Configuration")
          << "AllHitToTracksterAssociatorsProducer: 'layerClustersByCollection' has "
          << layerClustersByCollection.size() << " entries but 'tracksterCollections' has "
          << tracksterCollections.size() << ". They must match one-to-one.";
    }
    for (const auto& tag : layerClustersByCollection) {
      layerClustersByCollectionTokens_.emplace_back(consumes<std::vector<reco::CaloCluster>>(tag));
    }
  }

  for (const auto& tracksterToken : tracksterCollectionTokens_) {
    produces<ticl::AssociationMap<ticl::mapWithFraction>>("hitTo" + tracksterToken.first);
    produces<ticl::AssociationMap<ticl::mapWithFraction>>(tracksterToken.first + "ToHit");
  }
}

template <typename HIT>
void AllHitToTracksterAssociatorsProducerT<HIT>::produce(edm::StreamID,
                                                         edm::Event& iEvent,
                                                         const edm::EventSetup&) const {
  using namespace edm;

  Handle<std::vector<reco::CaloCluster>> layer_clusters;
  iEvent.getByToken(layerClustersToken_, layer_clusters);

  if (!layer_clusters.isValid()) {
    edm::LogWarning("AllHitToTracksterAssociatorsProducer") << "Missing LayerCluster collection.";
    for (const auto& tracksterToken : tracksterCollectionTokens_) {
      iEvent.put(std::make_unique<ticl::AssociationMap<ticl::mapWithFraction>>(), "hitTo" + tracksterToken.first);
      iEvent.put(std::make_unique<ticl::AssociationMap<ticl::mapWithFraction>>(), tracksterToken.first + "ToHit");
    }
    return;
  }

  Handle<std::unordered_map<DetId, const unsigned int>> hitMap;
  iEvent.getByToken(hitMapToken_, hitMap);

  if (!iEvent.getHandle(hitsToken_)) {
    edm::LogWarning("AllHitToTracksterAssociatorsProducer") << "Missing edm::RefProdVector<RecHitCollection>.";
    for (const auto& tracksterToken : tracksterCollectionTokens_) {
      iEvent.put(std::make_unique<ticl::AssociationMap<ticl::mapWithFraction>>(), "hitTo" + tracksterToken.first);
      iEvent.put(std::make_unique<ticl::AssociationMap<ticl::mapWithFraction>>(), tracksterToken.first + "ToHit");
    }
    return;
  }

  // Protection against missing RecHitCollection
  const auto hits = iEvent.get(hitsToken_);
  for (std::size_t index = 0; const auto& hitCollection : hits) {
    if (hitCollection->empty()) {
      LogDebug("AllHitToTracksterAssociatorsProducer") << "RecHitCollections #" << index << " is empty.";
    }
    index++;
  }

  edm::MultiSpan<HIT> rechitSpan(hits);
  // Check if rechitSpan is empty
  if (rechitSpan.size() == 0) {
    LogDebug("HitToSimClusterCaloParticleAssociatorProducer")
        << "Only empty RecHitCollections found. Association maps will be empty.";
    for (const auto& tracksterToken : tracksterCollectionTokens_) {
      iEvent.put(std::make_unique<ticl::AssociationMap<ticl::mapWithFraction>>(), "hitTo" + tracksterToken.first);
      iEvent.put(std::make_unique<ticl::AssociationMap<ticl::mapWithFraction>>(), tracksterToken.first + "ToHit");
    }
    return;
  }

  for (unsigned int ic = 0; ic < tracksterCollectionTokens_.size(); ++ic) {
    const auto& tracksterToken = tracksterCollectionTokens_[ic];
    Handle<std::vector<ticl::Trackster>> tracksters;
    iEvent.getByToken(tracksterToken.second, tracksters);

    if (!tracksters.isValid()) {
      LogDebug("AllHitToTracksterAssociatorsProducer")
          << "Missing Tracksters for collection " << tracksterToken.first << ". Association maps will be empty.";
      iEvent.put(std::make_unique<ticl::AssociationMap<ticl::mapWithFraction>>(), "hitTo" + tracksterToken.first);
      iEvent.put(std::make_unique<ticl::AssociationMap<ticl::mapWithFraction>>(), tracksterToken.first + "ToHit");
      continue;
    }

    // Resolve this collection's vertices against its own layer clusters when a
    // per-collection collection was configured, otherwise the shared one.
    const std::vector<reco::CaloCluster>* layerClustersForCollection = layer_clusters.product();
    if (!layerClustersByCollectionTokens_.empty()) {
      Handle<std::vector<reco::CaloCluster>> lcByCollection;
      iEvent.getByToken(layerClustersByCollectionTokens_[ic], lcByCollection);
      if (!lcByCollection.isValid()) {
        edm::LogWarning("AllHitToTracksterAssociatorsProducer")
            << "Missing per-collection LayerClusters for " << tracksterToken.first << ".";
        iEvent.put(std::make_unique<ticl::AssociationMap<ticl::mapWithFraction>>(), "hitTo" + tracksterToken.first);
        iEvent.put(std::make_unique<ticl::AssociationMap<ticl::mapWithFraction>>(), tracksterToken.first + "ToHit");
        continue;
      }
      layerClustersForCollection = lcByCollection.product();
    }

    auto hitToTracksterMap = std::make_unique<ticl::AssociationMap<ticl::mapWithFraction>>(rechitSpan.size());
    auto tracksterToHitMap = std::make_unique<ticl::AssociationMap<ticl::mapWithFraction>>(tracksters->size());

    for (unsigned int tracksterId = 0; tracksterId < tracksters->size(); ++tracksterId) {
      const auto& trackster = (*tracksters)[tracksterId];
      for (unsigned int j = 0; j < trackster.vertices().size(); ++j) {
        const auto& lc = (*layerClustersForCollection)[trackster.vertices()[j]];
        float invMultiplicity = 1.0f / trackster.vertex_multiplicity()[j];

        for (const auto& hitAndFraction : lc.hitsAndFractions()) {
          auto hitMapIter = hitMap->find(hitAndFraction.first);
          if (hitMapIter != hitMap->end()) {
            unsigned int rechitIndex = hitMapIter->second;
            float fraction = hitAndFraction.second * invMultiplicity;
            hitToTracksterMap->insert(rechitIndex, tracksterId, fraction);
            tracksterToHitMap->insert(tracksterId, rechitIndex, fraction);
          }
        }
      }
    }

    iEvent.put(std::move(hitToTracksterMap), "hitTo" + tracksterToken.first);
    iEvent.put(std::move(tracksterToHitMap), tracksterToken.first + "ToHit");
  }
}

template <typename HIT>
void AllHitToTracksterAssociatorsProducerT<HIT>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalMergeLayerClusters"));
  // Optional per-trackster-collection layer clusters (parallel to
  // tracksterCollections). Empty by default => use the single 'layerClusters'
  // for every collection.
  desc.add<std::vector<edm::InputTag>>("layerClustersByCollection", {});
  if constexpr (std::is_same_v<HIT, HGCRecHit>) {
    desc.add<std::vector<edm::InputTag>>("tracksterCollections",
                                         {edm::InputTag("ticlTrackstersCLUE3DHigh"),
                                          edm::InputTag("ticlTrackstersLinks"),
                                          edm::InputTag("ticlCandidate")});
    desc.add<edm::InputTag>("hitMapTag", edm::InputTag("recHitMapProducer", "hgcalRecHitMap"));
    desc.add<edm::InputTag>("hits", edm::InputTag("recHitMapProducer", "RefProdVectorHGCRecHitCollection"));
    descriptions.add("AllHitToTracksterAssociatorsProducer", desc);
  } else if constexpr (std::is_same_v<HIT, reco::PFRecHit>) {
    desc.add<std::vector<edm::InputTag>>("tracksterCollections", {edm::InputTag("ticlTrackstersCLUE3DBarrel")});
    desc.add<edm::InputTag>("hitMapTag", edm::InputTag("recHitMapProducer", "barrelRecHitMap"));
    desc.add<edm::InputTag>("hits", edm::InputTag("recHitMapProducer", "RefProdVectorPFRecHitCollection"));
    descriptions.add("AllHitToBarrelTracksterAssociatorsProducer", desc);
  }
}

template class AllHitToTracksterAssociatorsProducerT<HGCRecHit>;
template class AllHitToTracksterAssociatorsProducerT<reco::PFRecHit>;

// Define this as a plug-in
using AllHitToTracksterAssociatorsProducer = AllHitToTracksterAssociatorsProducerT<HGCRecHit>;
DEFINE_FWK_MODULE(AllHitToTracksterAssociatorsProducer);
using AllHitToBarrelTracksterAssociatorsProducer = AllHitToTracksterAssociatorsProducerT<reco::PFRecHit>;
DEFINE_FWK_MODULE(AllHitToBarrelTracksterAssociatorsProducer);
