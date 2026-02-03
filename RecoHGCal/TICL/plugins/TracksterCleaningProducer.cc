// Author: Riley Clark - riley.coltrane.clark@cern.ch
// Date: 10/2025

#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

#include "RecoHGCal/TICL/interface/TracksterCleaningAlgoBase.h"
#include "RecoHGCal/TICL/plugins/TracksterCleaningPluginFactory.h"

// Only needed for fillDescriptions (to populate parameters)
#include "RecoHGCal/TICL/plugins/TracksterCleaningByBeta.h"

using namespace ticl;

class TracksterCleaningProducer : public edm::stream::EDProducer<> {
public:
  explicit TracksterCleaningProducer(const edm::ParameterSet& ps);
  ~TracksterCleaningProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void produce(edm::Event& ev, const edm::EventSetup& es) override;

private:
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;
  edm::EDGetTokenT<std::vector<Trackster>> linked_token_;
  edm::EDGetTokenT<std::vector<Trackster>> clue3d_token_;

  edm::EDGetTokenT<std::vector<std::vector<unsigned int>>> map_token_;

  std::unique_ptr<TracksterCleaningAlgoBase> cleaningAlgo_;
  int algoVerbosity_{0};

  // cached for a clean one-time print
  edm::InputTag linkedTag_, clue3dTag_, layerClusterTag_, mapTag_;
};

TracksterCleaningProducer::TracksterCleaningProducer(const edm::ParameterSet& ps)
    : linkedTag_(ps.getParameter<edm::InputTag>("linkedTracksters")),
      clue3dTag_(ps.getParameter<edm::InputTag>("clue3DTracksters")),
      layerClusterTag_(ps.getParameter<edm::InputTag>("layer_clusters")),
      mapTag_(ps.getParameter<edm::InputTag>("clue3DInLinkedIndices")) {
  linked_token_   = consumes<std::vector<Trackster>>(linkedTag_);
  clue3d_token_   = consumes<std::vector<Trackster>>(clue3dTag_);
  clusters_token_ = consumes<std::vector<reco::CaloCluster>>(layerClusterTag_);
  map_token_      = consumes<std::vector<std::vector<unsigned int>>>(mapTag_);

  algoVerbosity_ = ps.getParameter<int>("algo_verbosity");

  const auto& cleanerPSet = ps.getParameter<edm::ParameterSet>("cleaner");
  const auto pluginName   = cleanerPSet.getParameter<std::string>("type");
  cleaningAlgo_ = std::unique_ptr<TracksterCleaningAlgoBase>(
      TracksterCleaningPluginFactory::get()->create(pluginName, cleanerPSet, consumesCollector()));

  produces<std::vector<Trackster>>();
  produces<std::vector<std::vector<unsigned int>>>();
  produces<std::vector<std::vector<unsigned int>>>("linkedTracksterIdToInputTracksterId");
}

void TracksterCleaningProducer::produce(edm::Event& ev, const edm::EventSetup& es) {
  static unsigned long long evtCount = 0;
  ++evtCount;

  const std::string label = moduleDescription().moduleLabel();

  auto const& linked        = ev.get(linked_token_);
  auto const& clue3d        = ev.get(clue3d_token_);
  auto const& layerClusters = ev.get(clusters_token_);
  auto const& mapIn         = ev.get(map_token_);

  // Skeletons fixed => must match exactly
  if (mapIn.size() != linked.size()) {
    throw cms::Exception("TracksterCleaningProducer")
        << "Size mismatch: linked.size()=" << linked.size()
        << " map.size()=" << mapIn.size()
        << " (map tag=" << mapTag_.encode() << ", linked tag=" << linkedTag_.encode() << ")";
  }

  // Guard: map indices must be valid clue3d indices
  for (unsigned int L = 0; L < mapIn.size(); ++L) {
    for (auto idx : mapIn[L]) {
      if (idx >= clue3d.size()) {
        throw cms::Exception("TracksterCleaningProducer")
            << "Invalid map index: map[" << L << "] has idx=" << idx
            << " clue3d.size()=" << clue3d.size();
      }
    }
  }

  // Guard: layer cluster indices referenced by clue3d must be valid
  for (unsigned int i = 0; i < clue3d.size(); ++i) {
    for (auto lcIdx : clue3d[i].vertices()) {
      if (lcIdx >= layerClusters.size()) {
        throw cms::Exception("TracksterCleaningProducer")
            << "Invalid layerCluster index: clue3d[" << i << "] has lcIdx=" << lcIdx
            << " layerClusters.size()=" << layerClusters.size();
      }
    }
  }

  auto outLinked = std::make_unique<std::vector<Trackster>>();
  auto outMap    = std::make_unique<std::vector<std::vector<unsigned int>>>();

  TracksterCleaningAlgoBase::Inputs in(ev, es, linked, clue3d, layerClusters, mapIn);
  cleaningAlgo_->cleanTracksters(in, *outLinked, *outMap);

  if (evtCount == 1 && algoVerbosity_ > 0) {
    std::cout << "[TICL-CLEAN][" << label << "] outputs:"
              << " outLinked=" << outLinked->size()
              << " outMap=" << outMap->size()
              << "\n";
  }

  // Keep identity links product (unchanged)
  auto outLinksDefault = std::make_unique<std::vector<std::vector<unsigned int>>>();
  outLinksDefault->resize(outLinked->size());
  for (unsigned int i = 0; i < outLinked->size(); ++i) {
    (*outLinksDefault)[i] = {i};
  }

  ev.put(std::move(outLinked));
  ev.put(std::move(outLinksDefault));
  ev.put(std::move(outMap), "linkedTracksterIdToInputTracksterId");
}

void TracksterCleaningProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("linkedTracksters", edm::InputTag("tracksterLinksProducer"));
  desc.add<edm::InputTag>("clue3DTracksters", edm::InputTag("ticlTrackstersCLUE3DHigh"));
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalMergeLayerClusters"));

  // IMPORTANT: tracked InputTag
  desc.add<edm::InputTag>("clue3DInLinkedIndices",
                          edm::InputTag("tracksterLinksProducer", "linkedTracksterIdToInputTracksterId"));

  desc.add<int>("algo_verbosity", 0);

  edm::ParameterSetDescription cleanerDesc;
  cleanerDesc.add<std::string>("type", "Beta");



  TracksterCleaningByBeta::fillPSetDescription(cleanerDesc);
  desc.add<edm::ParameterSetDescription>("cleaner", cleanerDesc);

  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(TracksterCleaningProducer);
