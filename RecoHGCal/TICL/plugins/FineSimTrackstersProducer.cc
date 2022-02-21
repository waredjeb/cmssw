// Author: Felice Pantaleo,Wahid Redjeb - felice.pantaleo@cern.ch,wahid.redjeb@cern.ch
// Date: 02/2022

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
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "DataFormats/HGCalReco/interface/TICLSeedingRegion.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "RecoHGCal/TICL/plugins/PatternRecognitionPluginFactory.h"
#include "PatternRecognitionbyCA.h"
#include "PatternRecognitionbyMultiClusters.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

using namespace ticl;

class FineSimTrackstersProducer : public edm::stream::EDProducer<edm::GlobalCache<TrackstersCache>> {
public:
  explicit FineSimTrackstersProducer(const edm::ParameterSet&, const TrackstersCache*);
  ~FineSimTrackstersProducer() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event&, const edm::EventSetup&) override;

  // static methods for handling the global cache
  static std::unique_ptr<TrackstersCache> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(TrackstersCache*);
  void computeSingleTracksterMask(std::vector<float>& result_mask, const Trackster& tracksters, const std::vector<reco::CaloCluster>& lcs){
    auto N = tracksters.vertices().size();
    for(size_t i_lc = 0; i_lc < N; i_lc++){
    result_mask[tracksters.vertices(i_lc)] = 1; // make the lcs available
    // std::cout << "SimTrackster index " << tracksters.vertices(i_lc) << std::endl;
    }
}

private:
  std::string detector_;
  bool doNose_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTiles>> myAlgo_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>> myAlgoHFNose_;
  const edm::EDGetTokenT<std::vector<Trackster>> simtrackster_token_;
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;
  const edm::EDGetTokenT<std::vector<float>> filtered_layerclusters_mask_token_;
  const edm::EDGetTokenT<std::vector<float>> original_layerclusters_mask_token_;
  const edm::EDGetTokenT<std::vector<TrackingParticle>> trkparticles_token_;
  const edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> clustersTime_token_;
  edm::EDGetTokenT<TICLLayerTiles> layer_clusters_tiles_token_;
  edm::EDGetTokenT<TICLLayerTilesHFNose> layer_clusters_tiles_hfnose_token_;
  const edm::EDGetTokenT<std::vector<TICLSeedingRegion>> seeding_regions_token_;
  const std::string itername_;
  ticl::Trackster::IterationIndex iterIndex_ = ticl::Trackster::IterationIndex(0);
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geom_token_;
  hgcal::RecHitTools rhtools_;
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
      simtrackster_token_(consumes<std::vector<Trackster>>(ps.getParameter<edm::InputTag>("simTracksters"))),
      clusters_token_(consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_clusters"))),
      filtered_layerclusters_mask_token_(consumes<std::vector<float>>(ps.getParameter<edm::InputTag>("filtered_mask"))),
      original_layerclusters_mask_token_(consumes<std::vector<float>>(ps.getParameter<edm::InputTag>("original_mask"))),
      trkparticles_token_(consumes(ps.getParameter<edm::InputTag>("trkparticles"))),
      clustersTime_token_(
          consumes<edm::ValueMap<std::pair<float, float>>>(ps.getParameter<edm::InputTag>("time_layerclusters"))),
      seeding_regions_token_(
          consumes<std::vector<TICLSeedingRegion>>(ps.getParameter<edm::InputTag>("seeding_regions"))),
      itername_(ps.getParameter<std::string>("itername")),
      geom_token_(esConsumes()) {
  auto plugin = ps.getParameter<std::string>("patternRecognitionBy");
  auto pluginPSet = ps.getParameter<edm::ParameterSet>("pluginPatternRecognitionBy" + plugin);
  if (doNose_) {
    myAlgoHFNose_ = PatternRecognitionHFNoseFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionBy"), pluginPSet, cache, consumesCollector());
  } else {
    myAlgo_ = PatternRecognitionFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionBy"), pluginPSet, cache, consumesCollector());
  }

  if (itername_ == "TrkEM")
    iterIndex_ = ticl::Trackster::TRKEM;
  else if (itername_ == "EM")
    iterIndex_ = ticl::Trackster::EM;
  else if (itername_ == "Trk")
    iterIndex_ = ticl::Trackster::TRKHAD;
  else if (itername_ == "HAD")
    iterIndex_ = ticl::Trackster::HAD;
  else if (itername_ == "MIP")
    iterIndex_ = ticl::Trackster::MIP;

  produces<std::vector<Trackster>>("fine");
  produces<std::vector<int>>("tracksterSeeds");
  produces<std::vector<std::vector<int>>>("tracksterSeedsDoublets");
  produces<std::vector<float>>("fine");  // Mask to be applied at the next iteration
}

void FineSimTrackstersProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // hgcalMultiClusters
  edm::ParameterSetDescription desc;
  desc.add<std::string>("detector", "HGCAL");
  desc.add<edm::InputTag>("simTracksters", edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<edm::InputTag>("filtered_mask", edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>("time_layerclusters", edm::InputTag("hgcalLayerClusters", "timeLayerCluster"));
  desc.add<edm::InputTag>("original_mask", edm::InputTag("hgcalLayerClusters", "InitialLayerClustersMask"));
  desc.add<edm::InputTag>("trkparticles", edm::InputTag("prunedTrackingParticles"));
  desc.add<edm::InputTag>("seeding_regions", edm::InputTag("ticlSeedingGlobal"));
  desc.add<std::string>("patternRecognitionBy", "CLUE3D");
  desc.add<std::string>("eid_graph_path", "RecoHGCal/TICL/data/tf_models/energy_id_v0.pb");
  desc.add<std::string>("itername", "unknown");

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

  descriptions.add("fineSimTrackstersProducer", desc);
}

void FineSimTrackstersProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
  std::cout << " --- BEGIN FINE SIMTRACKSTER PRODUCER ----- " << std::endl;
  auto result = std::make_unique<std::vector<Trackster>>();
  auto tracksterSeeds = std::make_unique<std::vector<int>>();
  auto tracksterSeedsDoublets = std::make_unique<std::vector<std::vector<int>>>();
  auto layer_clusters_tiles = std::make_unique<TICLLayerTiles>();
  auto layer_clusters_hfnose_tiles = std::make_unique<TICLLayerTilesHFNose>();
  // auto layer_clusters_tiles = std::make_unique<TICLLayerTiles>();
  // auto layer_clusters_hfnose_tiles = std::make_unique<TICLLayerTilesHFNose>();

  auto output_mask = std::make_unique<std::vector<float>>();
  const std::vector<Trackster>& simTracksters = evt.get(simtrackster_token_);
  const std::vector<float>& original_layerclusters_mask = evt.get(original_layerclusters_mask_token_);
  const auto& layerClusters = evt.get(clusters_token_);
  const auto& inputClusterMask = evt.get(filtered_layerclusters_mask_token_);
  const auto& layerClustersTimes = evt.get(clustersTime_token_);
  const auto& seeding_regions = evt.get(seeding_regions_token_);
  const auto& trkparticles = evt.get(trkparticles_token_);

  std::vector<float> fine_input_cluster_mask(layerClusters.size(), 0);
  std::vector<float> cluters_mask(layerClusters.size(), 0.);
  // fine_input_cluster_mask.reserve(layerClusters.size());
  const auto& geom = es.getData(geom_token_);
  rhtools_.setGeometry(geom);
  double counter = 0.;
  int lcId = 0;

  for (size_t i = 0; i < inputClusterMask.size(); i++) {
    if (inputClusterMask[i] == 0) {
      fine_input_cluster_mask[i] = 1;
      counter += 1;
      auto lc = layerClusters[i];
      const auto firstHitDetId = lc.hitsAndFractions()[0].first;
      int layer = rhtools_.getLayerWithOffset(firstHitDetId) +
                  rhtools_.lastLayer(doNose_) * ((rhtools_.zside(firstHitDetId) + 1) >> 1) - 1;
      assert(layer >= 0);

      if (doNose_)
        layer_clusters_hfnose_tiles->fill(layer, lc.eta(), lc.phi(), lcId);
      else
        layer_clusters_tiles->fill(layer, lc.eta(), lc.phi(), lcId);
      lcId++;
    } else {
      fine_input_cluster_mask[i] = 0;
    }
  }

    // std::vector<float> fine_input_cluster_mask;
    auto tmp_result = std::make_unique<std::vector<Trackster>>();
    // fine_input_cluster_mask.resize(layerClusters.size(), 0);
    // computeSingleTracksterMask(fine_input_cluster_mask, simTracksters[i_st], layerClusters);
    auto count_av_lcs = std::count_if(fine_input_cluster_mask.begin(), fine_input_cluster_mask.end(), [](int i ){return i == 1;});
    auto count_nonav_lcs = std::count_if(fine_input_cluster_mask.begin(), fine_input_cluster_mask.end(), [](int i ){return i == 0;});
    std::cout << "Count av lcs " << count_av_lcs << std::endl;
    std::cout << "Count non av lcs " << count_nonav_lcs << std::endl;

    std::unordered_map<int, std::vector<int>> seedToTrackstersAssociation;
    // if it's regional iteration and there are seeding regions
    if (!seeding_regions.empty() and seeding_regions[0].index != -1) {
      auto numberOfSeedingRegions = seeding_regions.size();
      for (unsigned int i = 0; i < numberOfSeedingRegions; ++i) {
        seedToTrackstersAssociation.emplace(seeding_regions[i].index, 0);
      }
    }
    if (doNose_) {
      const typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Inputs inputHFNose(evt,
                                                                                          es,
                                                                                          layerClusters,
                                                                                          fine_input_cluster_mask,
                                                                                          layerClustersTimes,
                                                                                          *layer_clusters_hfnose_tiles,
                                                                                          seeding_regions);

      typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Outputs output(
          *tmp_result, *tracksterSeeds, *tracksterSeedsDoublets);
      myAlgoHFNose_->makeTracksters(inputHFNose, output, seedToTrackstersAssociation);

    } else {
      const typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Inputs input(
          evt, es, layerClusters, fine_input_cluster_mask, layerClustersTimes, *layer_clusters_tiles, seeding_regions);

      typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Outputs output(
          *tmp_result, *tracksterSeeds, *tracksterSeedsDoublets);
      myAlgo_->makeTracksters(input, output, seedToTrackstersAssociation);
    }
  for(auto& tr : *tmp_result){
    std::cout << " RESULTS LCs size " << tr.vertices().size() << std::endl;
  }
  result->insert(result->end(), tmp_result->begin(), tmp_result->end());
  // }
  // Now update the global mask and put it into the event
  output_mask->reserve(original_layerclusters_mask.size());
  // Copy over the previous state
  std::copy(
      std::begin(original_layerclusters_mask), std::end(original_layerclusters_mask), std::back_inserter(*output_mask));
  
  std::cout << " SimTracksters Size " << simTracksters.size() <<  " FineSimTrackster Size " << result->size()  << std::endl;
  evt.put(std::move(result), "fine");
  evt.put(std::move(output_mask), "fine");
  // evt.put(std::move(result_fromCP), "fineFromCPs");
  // evt.put(std::move(output_mask_fromCP), "fineFromCPs");
  // evt.put(std::move(cpToSc_SimTrackstersMap), "fine");
  std::cout << " --- END FINE SIMTRACKSTER PRODUCER ----- " << std::endl;
}
