// Author: Felice Pantaleo, Wahid Redjeb - felice.pantaleo@cern.ch, wahid.redjeb@cern.ch
// Date: 02/2022

// user include files

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
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "SimDataFormats/Associations/interface/LayerClusterToSimClusterAssociator.h"
#include "SimDataFormats/Associations/interface/LayerClusterToCaloParticleAssociator.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "DataFormats/HGCalReco/interface/TICLSeedingRegion.h"

#include "RecoHGCal/TICL/plugins/PatternRecognitionPluginFactory.h"
#include "PatternRecognitionbyCA.h"

#include "RecoHGCal/TICL/interface/commons.h"

#include "TrackstersPCA.h"
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>

using namespace ticl;

class FineSimTrackstersProducer : public edm::stream::EDProducer<edm::GlobalCache<TrackstersCache>> {
public:
  explicit FineSimTrackstersProducer(const edm::ParameterSet&, const TrackstersCache*);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event&, const edm::EventSetup&) override;

  // static methods for handling the global cache
  static std::unique_ptr<TrackstersCache> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(TrackstersCache*);

private:
  std::string detector_;
  const bool doNose_ = false;
  const edm::EDGetTokenT<std::vector<Trackster>> simTracksters_token_;
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;
  const edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> clustersTime_token_;
  const edm::EDGetTokenT<std::vector<float>> filtered_layerclusters_mask_token_;

  const edm::EDGetTokenT<std::vector<SimCluster>> simclusters_token_;
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloparticles_token_;
  const edm::EDGetTokenT<std::vector<TrackingParticle>> trkparticles_token_;
  const edm::EDGetTokenT<std::vector<TICLSeedingRegion>> seeding_regions_token_;

  const edm::EDGetTokenT<hgcal::SimToRecoCollectionWithSimClusters> associatorMapSimClusterToReco_token_;
  const edm::EDGetTokenT<hgcal::SimToRecoCollection> associatorMapCaloParticleToReco_token_;
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geom_token_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTiles>> myAlgo_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>> myAlgoHFNose_;
  hgcal::RecHitTools rhtools_;

  const double fractionCut_;
};
DEFINE_FWK_MODULE(FineSimTrackstersProducer);

std::unique_ptr<TrackstersCache> FineSimTrackstersProducer::initializeGlobalCache(const edm::ParameterSet& params) {
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

void FineSimTrackstersProducer::globalEndJob(TrackstersCache* cache) {
  delete cache->eidGraphDef;
  cache->eidGraphDef = nullptr;
}

FineSimTrackstersProducer::FineSimTrackstersProducer(const edm::ParameterSet& ps, const TrackstersCache* cache)
    : detector_(ps.getParameter<std::string>("detector")),
      doNose_(detector_ == "HFNose"),
      simTracksters_token_(consumes(ps.getParameter<edm::InputTag>("simTracksters"))),
      clusters_token_(consumes(ps.getParameter<edm::InputTag>("layer_clusters"))),
      clustersTime_token_(consumes(ps.getParameter<edm::InputTag>("time_layerclusters"))),
      filtered_layerclusters_mask_token_(consumes(ps.getParameter<edm::InputTag>("filtered_mask"))), //consume the filtered one to run PR on SimTracksters' LCs
      simclusters_token_(consumes(ps.getParameter<edm::InputTag>("simclusters"))),
      caloparticles_token_(consumes(ps.getParameter<edm::InputTag>("caloparticles"))),
      trkparticles_token_(consumes(ps.getParameter<edm::InputTag>("trkparticles"))),
      seeding_regions_token_(
          consumes<std::vector<TICLSeedingRegion>>(ps.getParameter<edm::InputTag>("seeding_regions"))),
      associatorMapSimClusterToReco_token_(
          consumes(ps.getParameter<edm::InputTag>("layerClusterSimClusterAssociator"))),
      associatorMapCaloParticleToReco_token_(
          consumes(ps.getParameter<edm::InputTag>("layerClusterCaloParticleAssociator"))),
      geom_token_(esConsumes()),
      fractionCut_(ps.getParameter<double>("fractionCut")) {
  auto plugin = ps.getParameter<std::string>("patternRecognitionBy");
  auto pluginPSet = ps.getParameter<edm::ParameterSet>("pluginPatternRecognitionBy" + plugin);
  if (doNose_) {
    myAlgoHFNose_ = PatternRecognitionHFNoseFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionBy"), pluginPSet, cache, consumesCollector());
  } else {
    myAlgo_ = PatternRecognitionFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionBy"), pluginPSet, cache, consumesCollector());
  }
  produces<TracksterCollection>("fine");
  produces<std::vector<float>>("fine");
  produces<TracksterCollection>("fineFromCPs");
  produces<std::vector<float>>("fineFromCPs");
  produces<std::map<uint, std::vector<uint>>>("fine");
}

void FineSimTrackstersProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // hgcalMultiClusters
  edm::ParameterSetDescription desc;
  desc.add<std::string>("detector", "HGCAL");
  desc.add<edm::InputTag>("simTracksters", edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<edm::InputTag>("time_layerclusters", edm::InputTag("hgcalLayerClusters", "timeLayerCluster"));
  desc.add<edm::InputTag>("filtered_mask", edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>("simclusters", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("caloparticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("trkparticles", edm::InputTag("prunedTrackingParticles"));
  desc.add<edm::InputTag>("seeding_regions", edm::InputTag("ticlSeedingGlobal"));
  desc.add<std::string>("patternRecognitionBy", "CA");
  desc.add<edm::InputTag>("layerClusterSimClusterAssociator",
                          edm::InputTag("layerClusterSimClusterAssociationProducer"));
  desc.add<edm::InputTag>("layerClusterCaloParticleAssociator",
                          edm::InputTag("layerClusterCaloParticleAssociationProducer"));
  desc.add<double>("fractionCut", 0.);
  desc.add<std::string>("eid_graph_path", "RecoHGCal/TICL/data/tf_models/energy_id_v0.pb");

  // CA Plugin
  edm::ParameterSetDescription pluginDesc;
  pluginDesc.addNode(edm::PluginDescription<PatternRecognitionFactory>("type", "CA", true));
  desc.add<edm::ParameterSetDescription>("pluginPatternRecognitionByCA", pluginDesc);
  //
  // CLUE3D Plugin
  edm::ParameterSetDescription pluginDescClue3D;
  pluginDescClue3D.addNode(edm::PluginDescription<PatternRecognitionFactory>("type", "CLUE3D", true));
  desc.add<edm::ParameterSetDescription>("pluginPatternRecognitionByCLUE3D", pluginDescClue3D);

  // FastJet Plugin
  edm::ParameterSetDescription pluginDescFastJet;
  pluginDescFastJet.addNode(edm::PluginDescription<PatternRecognitionFactory>("type", "FastJet", true));
  desc.add<edm::ParameterSetDescription>("pluginPatternRecognitionByFastJet", pluginDescFastJet);

  descriptions.addWithDefaultLabel(desc);
}

void FineSimTrackstersProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
  auto tiles = std::make_unique<TICLLayerTiles>();
  auto tilesHFNose = std::make_unique<TICLLayerTilesHFNose>();
  auto result = std::make_unique<TracksterCollection>();
  auto output_mask = std::make_unique<std::vector<float>>();
  auto result_fromCP = std::make_unique<TracksterCollection>();
  auto output_mask_fromCP = std::make_unique<std::vector<float>>();
  auto cpToSc_SimTrackstersMap = std::make_unique<std::map<uint, std::vector<uint>>>();
  auto tracksterSeeds = std::make_unique<std::vector<int>>();
  auto tracksterSeedsDoublets = std::make_unique<std::vector<std::vector<int>>>();

  const auto& simTracksters = evt.get(simTracksters_token_);
  const auto& layerClusters = evt.get(clusters_token_);
  const auto& layerClustersTimes = evt.get(clustersTime_token_);
  const auto& inputClusterMask = evt.get(filtered_layerclusters_mask_token_);
  output_mask->resize(layerClusters.size(), 1.f);
  output_mask_fromCP->resize(layerClusters.size(), 1.f);

  const auto& simclusters = evt.get(simclusters_token_);
  const auto& caloparticles = evt.get(caloparticles_token_);
  const auto& trkparticles = evt.get(trkparticles_token_);

  const auto& simClustersToRecoColl = evt.get(associatorMapSimClusterToReco_token_);
  const auto& caloParticlesToRecoColl = evt.get(associatorMapCaloParticleToReco_token_);

  const auto& geom = es.getData(geom_token_);
  rhtools_.setGeometry(geom);
  const auto num_simclusters = simclusters.size();
  // result->reserve(num_simclusters);  // Conservative size, will call shrink_to_fit later
  const auto num_caloparticles = caloparticles.size();
  // result_fromCP->reserve(num_caloparticles);

  auto sim_tracksters_size = simTracksters.size();
  std::vector<TICLSeedingRegion> seeding_regions;
  for (size_t i = 0; i < sim_tracksters_size; i++) {
    auto sim_t = simTracksters[i];
    auto N_lcs = sim_t.vertices().size();
    for (size_t i_lc = 0; i_lc < N_lcs; i_lc++) {
      auto lc = layerClusters[sim_t.vertices(i_lc)];
    }
  }

  for (size_t i = 0; i < layerClusters.size(); i++) {
    if (inputClusterMask[i] > 0) {
      auto lc = layerClusters[i];
      const auto firstHitDetId = lc.hitsAndFractions()[0].first;
      int layer = rhtools_.getLayerWithOffset(firstHitDetId) +
                  rhtools_.lastLayer(doNose_) * ((rhtools_.zside(firstHitDetId) + 1) >> 1) - 1;
      assert(layer >= 0);

      if (doNose_)
        tilesHFNose->fill(layer, lc.eta(), lc.phi(), i);
      else
        tiles->fill(layer, lc.eta(), lc.phi(), i);
    }
  }

  std::unordered_map<int, std::vector<int>> seedToTrackstersAssociation;
  // if it's regional iteration and there are seeding regions
  // nothing done here, no seeding regions, just global iterations!
  if (!seeding_regions.empty() and seeding_regions[0].index != -1) {
    auto numberOfSeedingRegions = seeding_regions.size();
    for (unsigned int i = 0; i < numberOfSeedingRegions; ++i) {
      seedToTrackstersAssociation.emplace(seeding_regions[i].index, 0);
    }
  }

  if (doNose_) {
    // const auto& layer_clusters_hfnose_tiles = evt.get(layer_clusters_tiles_hfnose_token_);
    const typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Inputs inputHFNose(
        evt, es, layerClusters, inputClusterMask, layerClustersTimes, *(tilesHFNose.get()), seeding_regions);

    typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Outputs output(
        *result, *tracksterSeeds, *tracksterSeedsDoublets);
    myAlgoHFNose_->makeTracksters(inputHFNose, output, seedToTrackstersAssociation);

  } else {
    // const auto& layer_clusters_tiles = evt.get(layer_clusters_tiles_token_);
    const typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Inputs input(
        evt, es, layerClusters, inputClusterMask, layerClustersTimes, *(tiles.get()), seeding_regions);
    typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Outputs output(
        *result, *tracksterSeeds, *tracksterSeedsDoublets);
    myAlgo_->makeTracksters(input, output, seedToTrackstersAssociation);
  }

  evt.put(std::move(result), "fine");
  evt.put(std::move(output_mask), "fine");
  evt.put(std::move(result_fromCP), "fineFromCPs");
  evt.put(std::move(output_mask_fromCP), "fineFromCPs");
  evt.put(std::move(cpToSc_SimTrackstersMap), "fine");
}
