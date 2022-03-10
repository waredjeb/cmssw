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
#include "RecoHGCal/TICL/plugins/TrackstersPCA.h"
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
  void computeSingleTracksterMask(std::vector<float>& result_mask,
                                  std::vector<float>& test,
                                  const Trackster& trackster,
                                  const std::vector<reco::CaloCluster>& lcs) {
    auto N = trackster.vertices().size();
    for (size_t i_lc = 0; i_lc < N; i_lc++) {
      result_mask[trackster.vertices(i_lc)] =
          1. / trackster.vertex_multiplicity(i_lc);  // make the lcs available weighted by the vertex_multiplicity.
      test[trackster.vertices(i_lc)] = 1.;
    }
  }
  int findNearestSeed(const reco::CaloCluster& lc, std::vector<int>& tracksterSeeds, const std::vector<reco::CaloCluster>& lcs){
    double min = 1e9;
    int result = -1;
    auto distance_sq = [](reco::CaloCluster lc1,  reco::CaloCluster lc2){
      auto dx = lc2.x() - lc1.x();
      auto dy = lc2.y() - lc1.y();
      auto dz = lc2.z() - lc1.z();
      auto dx2 = dx*dx;
      auto dy2 = dy*dy;
      auto dz2 = dz*dz;

      return sqrt(dx2 + dy2 + dz2);
    };
    auto r = 0;
    for(auto & i_s : tracksterSeeds){
      auto distance = distance_sq(lc, lcs[i_s]);
      if(distance < min){
        min = distance;
        result = r;
      }
      r++;
    }
    return result;
  }

  // static methods for handling the global cache
  static std::unique_ptr<TrackstersCache> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(TrackstersCache*);

private:
  std::string detector_;
  bool doNose_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTiles>> myAlgo_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>> myAlgoHFNose_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTiles>> myAlgoMIP_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>> myAlgoHFNoseMIP_;
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
  auto pluginMIP = ps.getParameter<std::string>("patternRecognitionMIPBy");
  auto pluginPSet = ps.getParameter<edm::ParameterSet>("pluginPatternRecognitionBy" + plugin);
  auto pluginPSetMIP = ps.getParameter<edm::ParameterSet>("pluginPatternRecognitionMIPBy" + pluginMIP);
  if (doNose_) {
    myAlgoHFNose_ = PatternRecognitionHFNoseFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionBy"), pluginPSet, cache, consumesCollector());
    layer_clusters_tiles_hfnose_token_ =
        consumes<TICLLayerTilesHFNose>(ps.getParameter<edm::InputTag>("layer_clusters_hfnose_tiles"));
    myAlgoHFNoseMIP_ = PatternRecognitionHFNoseFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionMIPBy"), pluginPSetMIP, cache, consumesCollector());
  } else {
    myAlgo_ = PatternRecognitionFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionBy"), pluginPSet, cache, consumesCollector());
    layer_clusters_tiles_token_ = consumes<TICLLayerTiles>(ps.getParameter<edm::InputTag>("layer_clusters_tiles"));
    myAlgoMIP_ = PatternRecognitionFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionMIPBy"), pluginPSetMIP, cache, consumesCollector());
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
  produces<std::vector<int>>("fine");
  produces<std::vector<Trackster>>("fineMIP");
  produces<std::vector<int>>("fineMIP");
  produces<std::map<uint, std::vector<uint>>>("fine");
  // produces<std::vector<std::vector<int>>>("tracksterSeedsDoublets");
  produces<std::vector<float>>("fine");  //  Mask to be applied at the next iteration
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
  desc.add<edm::InputTag>("layer_clusters_tiles", edm::InputTag("ticlLayerTileProducer"));
  desc.add<edm::InputTag>("layer_clusters_hfnose_tiles", edm::InputTag("ticlLayerTileHFNose"));
  desc.add<std::string>("patternRecognitionBy", "CLUE3D");
  desc.add<std::string>("patternRecognitionMIPBy", "CA");
  desc.add<std::string>("eid_graph_path", "RecoHGCal/TICL/data/tf_models/energy_id_v0.pb");
  desc.add<std::string>("itername", "unknown");

  // CA Plugin
  edm::ParameterSetDescription pluginDesc;
  pluginDesc.addNode(edm::PluginDescription<PatternRecognitionFactory>("type", "CA", true));
  desc.add<edm::ParameterSetDescription>("pluginPatternRecognitionMIPByCA", pluginDesc);
  //CA Plugin for MIP
  // edm::ParameterSetDescription pluginDescCA;
  // pluginDescCA.addNode(edm::PluginDescription<PatternRecognitionFactory>("type", "CA", true));
  // desc.add<edm::ParameterSetDescription>("pluginPatternRecognitionByCA", pluginDescCA);
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
  std::cout << " BEGIN PRODUCER " << std::endl;
  auto result = std::make_unique<std::vector<Trackster>>();
  auto resultMIP = std::make_unique<std::vector<Trackster>>();
  auto tracksterSeedsMIP = std::make_unique<std::vector<int>>();
  auto tracksterSeeds = std::make_unique<std::vector<int>>();
  auto tracksterSeedsDoublets = std::make_unique<std::vector<std::vector<int>>>();
  auto layer_clusters_tiles = std::make_unique<TICLLayerTiles>();
  auto layer_clusters_hfnose_tiles = std::make_unique<TICLLayerTilesHFNose>();
  auto simTracksterToFineSimTracksters = std::make_unique<std::map<uint, std::vector<uint>>>();
  const std::vector<Trackster>& simTracksters = evt.get(simtrackster_token_);
  const auto& layerClusters = evt.get(clusters_token_);
  const auto& layerClustersTimes = evt.get(clustersTime_token_);
  const auto& seeding_regions = evt.get(seeding_regions_token_);
  const auto& geom = es.getData(geom_token_);
  rhtools_.setGeometry(geom);
  double counter = 0.;
  int lcId = 0;
  int t = 0;
  std::cout << "Last Layer EE z" << rhtools_.getPositionLayer(rhtools_.lastLayerEE(doNose_)).z() << std::endl;
  std::cout << "Last Layer FG z" << rhtools_.getPositionLayer(rhtools_.lastLayerFH()).z() << std::endl;

  auto index_tmp_finesimts = 0;

  for (size_t i_st = 0; i_st < simTracksters.size(); i_st++) {
    int tot_number_of_lcs_in_fineSimTracksters = 0;
    int tot_number_of_lcs_in_simTracksters = 0;
    tot_number_of_lcs_in_simTracksters += simTracksters[i_st].vertices().size();
    std::vector<uint> fine_sim_trackster_index;
    std::vector<float> fine_input_cluster_mask(layerClusters.size(), 0.);
    std::vector<float> test_fine_input_cluster_mask(layerClusters.size(), 0.);
    computeSingleTracksterMask(fine_input_cluster_mask, test_fine_input_cluster_mask, simTracksters[i_st], layerClusters);
    std::vector<float> output_mask(fine_input_cluster_mask);
    auto count_av = 0;
    auto count_av_2 = std::count_if(fine_input_cluster_mask.begin(), fine_input_cluster_mask.end(), [](float i ) {return i > 0.;});
    for(size_t f_i = 0; f_i < fine_input_cluster_mask.size(); ++f_i){
      if(fine_input_cluster_mask[f_i] > 0.){
        count_av += 1;
      }
    }

    if(count_av != static_cast<int>(simTracksters[i_st].vertices().size()) || count_av_2 != static_cast<int>(simTracksters[i_st].vertices().size())){
      std::cout << "DIFFERENCE " << count_av - static_cast<int>(simTracksters[i_st].vertices().size()) << std::endl;
      std::cout << "DIFFERENCE_2  " << count_av_2 - static_cast<int>(simTracksters[i_st].vertices().size()) << std::endl;
    
      for(size_t j = 0; j < simTracksters[i_st].vertices().size(); ++j){
        auto lc_id = simTracksters[i_st].vertices(j);
        if(!(fine_input_cluster_mask[lc_id] > 0.)){
          std::cout << " Mask " << fine_input_cluster_mask[lc_id] << " Vertex Multiplicity " <<  simTracksters[i_st].vertex_multiplicity(j) << std::endl;
        }
      }
    }
    auto tmp_result = std::make_unique<std::vector<Trackster>>();
    auto tmp_tracksterSeeds = std::make_unique<std::vector<int>>();
    auto tmp_resultMIP = std::make_unique<std::vector<Trackster>>();
    auto tmp_tracksterSeedsMIP = std::make_unique<std::vector<int>>();
    auto tmp_tracksterSeedsDoubletsMIP = std::make_unique<std::vector<std::vector<int>>>();
    auto tmp_tracksterSeedsDoublets = std::make_unique<std::vector<std::vector<int>>>();
    std::unordered_map<int, std::vector<int>> seedToTrackstersAssociation;
    // if it's regional iteration and there are seeding regions
    if (!seeding_regions.empty() and seeding_regions[0].index != -1) {
      auto numberOfSeedingRegions = seeding_regions.size();
      for (unsigned int i = 0; i < numberOfSeedingRegions; ++i) {
        seedToTrackstersAssociation.emplace(seeding_regions[i].index, 0);
      }
    }
    if (doNose_) {
      const auto& layer_clusters_hfnose_tiles = evt.get(layer_clusters_tiles_hfnose_token_);
      const typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Inputs inputHFNose(evt,
                                                                                           es,
                                                                                           layerClusters,
                                                                                           fine_input_cluster_mask,
                                                                                           layerClustersTimes,
                                                                                           layer_clusters_hfnose_tiles,
                                                                                           seeding_regions);

      typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Outputs output(
          *tmp_result, *tmp_tracksterSeeds, *tmp_tracksterSeedsDoublets);
      myAlgoHFNose_->makeTracksters(inputHFNose, output, seedToTrackstersAssociation);

    } else {
      const auto& layer_clusters_tiles = evt.get(layer_clusters_tiles_token_);
      const typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Inputs input(
          evt, es, layerClusters, fine_input_cluster_mask, layerClustersTimes, layer_clusters_tiles, seeding_regions);

      typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Outputs output(
          *tmp_result, *tmp_tracksterSeeds, *tmp_tracksterSeedsDoublets);
      myAlgo_->makeTracksters(input, output, seedToTrackstersAssociation);
    }

    for (auto& t_fst : *tmp_result) {
      for (size_t i_lc = 0; i_lc < t_fst.vertices().size(); ++i_lc) {
        auto lc_id = t_fst.vertices(i_lc);
      }
    }
    //if no seeds, just take the simtrackster.
    if (tmp_tracksterSeeds->empty()) {
      tmp_result->push_back(simTracksters[i_st]);
    }

    // LCs recovery
    for (auto& trackster : *tmp_result) {
      // Mask the used elements, accordingly
      for (auto const v : trackster.vertices()) {
        // TODO(rovere): for the moment we mask the layer cluster completely. In
        // the future, properly compute the fraction of usage.
        output_mask[v] = 0.;
      }
    }
  auto tot_still_av = 0;
  
  for(auto& x : output_mask){
    if(x > 0){
      tot_still_av +=1;
    }
  }
  std::cout << " LCs still available " << tot_still_av << std::endl;
  if (doNose_) {
      const auto& layer_clusters_hfnose_tiles_mip = evt.get(layer_clusters_tiles_hfnose_token_);
      const typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Inputs inputHFNoseMIP(evt,
                                                                                           es,
                                                                                           layerClusters,
                                                                                           output_mask,
                                                                                           layerClustersTimes,
                                                                                           layer_clusters_hfnose_tiles_mip,
                                                                                           seeding_regions);

      typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Outputs outputMIPNose(
          *tmp_resultMIP, *tmp_tracksterSeedsMIP, *tmp_tracksterSeedsDoubletsMIP);
      myAlgoHFNoseMIP_->makeTracksters(inputHFNoseMIP, outputMIPNose, seedToTrackstersAssociation);

    } else {
      const auto& layer_clusters_tiles_mip = evt.get(layer_clusters_tiles_token_);
      const typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Inputs inputMIP(
          evt, es, layerClusters, output_mask, layerClustersTimes, layer_clusters_tiles_mip, seeding_regions);

      typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Outputs outputMIP(
          *tmp_resultMIP, *tmp_tracksterSeedsMIP, *tmp_tracksterSeedsDoubletsMIP);
      myAlgoMIP_->makeTracksters(inputMIP, outputMIP, seedToTrackstersAssociation);
    }

    // for(size_t i = 0; i < simTracksters[i_st].vertices().size(); i++){
    //   auto i_lc_st = simTracksters[i_st].vertices(i);
    //   auto index_nearest_seed = -1;
    //   if(output_mask[i_lc_st]  > 0.){
    //     auto lc = layerClusters[i_lc_st];
    //     index_nearest_seed = findNearestSeed(lc, *tmp_tracksterSeeds, layerClusters);
    //     std::cout << "Index nearest seed " << index_nearest_seed << " Mask " << output_mask[i_lc_st] <<  std::endl;
    //     if(index_nearest_seed >= 0){
    //       (*tmp_result)[index_nearest_seed].vertices().push_back(i_lc_st);
    //     }
    //   }
    // }

    // (*simTracksterToFineSimTracksters)[i_st] = fine_sim_trackster_index;
    result->insert(result->end(), tmp_result->begin(), tmp_result->end());
    resultMIP->insert(resultMIP->end(), tmp_resultMIP->begin(), tmp_resultMIP->end());
    tracksterSeedsMIP->insert(tracksterSeedsMIP->end(), tmp_tracksterSeedsMIP->begin(), tmp_tracksterSeedsMIP->end());
    tracksterSeeds->insert(tracksterSeeds->end(), tmp_tracksterSeeds->begin(), tmp_tracksterSeeds->end());
  }  // end simtracksters loop

  ticl::assignPCAtoTracksters(*result,
                          layerClusters,
                          layerClustersTimes,
                          rhtools_.getPositionLayer(rhtools_.lastLayerEE(false), false).z());

  // run energy regression and ID
  // energyRegressionAndID(layerClusters, result);


  evt.put(std::move(result), "fine");
  evt.put(std::move(resultMIP), "fineMIP");
  evt.put(std::move(simTracksterToFineSimTracksters), "fine");
  evt.put(std::move(tracksterSeeds), "fine");
  evt.put(std::move(tracksterSeedsMIP), "fineMIP");
  std::cout << "END PRODUCER " << std::endl;
}
