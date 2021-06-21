// user include files
#include <vector>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"

#include "RecoHGCal/TICL/interface/TrackstersLinkingPluginFactory.h"
#include "TrackstersLinkingbyPCA.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

using namespace ticl;

class LinkedTrackstersProducer : public edm::stream::EDProducer<edm::GlobalCache<TrackstersCache>> {
public:
  explicit LinkedTrackstersProducer(const edm::ParameterSet&, const TrackstersCache*);
  ~LinkedTrackstersProducer() override {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event&, const edm::EventSetup&) override;

  // static methods for handling the global cache
  static std::unique_ptr<TrackstersCache> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(TrackstersCache*);

private:
  std::string detector_;
  std::unique_ptr<TrackstersLinkingAlgoBaseT> myAlgo_;
  const edm::EDGetTokenT<std::vector<SimCluster>> clusters_token_;
  const edm::EDGetTokenT<std::vector<Trackster>> tracksters_token_;
};
DEFINE_FWK_MODULE(LinkedTrackstersProducer);

std::unique_ptr<TrackstersCache> LinkedTrackstersProducer::initializeGlobalCache(const edm::ParameterSet& params) {
  // this method is supposed to create, initialize and return a TrackstersCache instance
  std::unique_ptr<TrackstersCache> cache = std::make_unique<TrackstersCache>(params);

  // load the graph def and save it
  std::string graphPath = params.getParameter<std::string>("eid_graph_path");
  if (!graphPath.empty()) {
    graphPath = edm::FileInPath(graphPath).fullPath();
    cache->eidGraphDef = tensorflow::loadGraphDef(graphPath);
  }

  return cache;
}

void LinkedTrackstersProducer::globalEndJob(TrackstersCache* cache) {
  delete cache->eidGraphDef;
  cache->eidGraphDef = nullptr;
}

LinkedTrackstersProducer::LinkedTrackstersProducer(const edm::ParameterSet& ps, const TrackstersCache* cache)
    : detector_(ps.getParameter<std::string>("detector")),
      clusters_token_(consumes<std::vector<SimCluster>>(ps.getParameter<edm::InputTag>("simCluster"))),
      tracksters_token_(consumes<std::vector<Trackster>>(ps.getParameter<edm::InputTag>("tracksters"))) {
  auto plugin = ps.getParameter<std::string>("trackstersLinkingBy");
  auto pluginPSet = ps.getParameter<edm::ParameterSet>("pluginTrackstersLinkingBy" + plugin);
  myAlgo_ =
      TrackstersLinkingFactory::get()->create(ps.getParameter<std::string>("trackstersLinkingBy"), pluginPSet, cache);
  produces<std::vector<Trackster>>();
  produces<std::vector<float>>();  // Mask to be applied at the next iteration
}

void LinkedTrackstersProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // hgcalMultiClusters
  edm::ParameterSetDescription desc;
  desc.add<std::string>("detector", "HGCAL");
  desc.add<edm::InputTag>("simCluster", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("tracksters", edm::InputTag("ticlSimTracksters"));
  desc.add<std::string>("trackstersLinkingBy", "PCA");
  desc.add<std::string>("eid_graph_path", "RecoHGCal/TICL/data/tf_models/energy_id_v0.pb");
  // PCA Plugin
  edm::ParameterSetDescription pluginDesc;
  pluginDesc.addNode(edm::PluginDescription<TrackstersLinkingFactory>("type", "PCA", true));
  desc.add<edm::ParameterSetDescription>("pluginTrackstersLinkingByPCA", pluginDesc);

  descriptions.add("linkedTrackstersProducer", desc);
}

void LinkedTrackstersProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
  auto result = std::make_unique<std::vector<float>>();
  const auto& layerClusters = evt.get(clusters_token_);
  const auto& tracksters = evt.get(tracksters_token_);
}