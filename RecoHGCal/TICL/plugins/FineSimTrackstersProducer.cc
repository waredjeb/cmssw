// Author: Felice Pantaleo,Marco Rovere - felice.pantaleo@cern.ch,marco.rovere@cern.ch
// Date: 09/2018

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

#include "RecoHGCal/TICL/plugins/PatternRecognitionPluginFactory.h"
#include "PatternRecognitionbyCA.h"
#include "PatternRecognitionbyMultiClusters.h"
#include "TrackstersPCA.h"

#include "PhysicsTools/TensorFlow/interface/TfGraphRecord.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "PhysicsTools/TensorFlow/interface/TfGraphDefWrapper.h"

using namespace ticl;

class FineSimFineSimTrackstersProducer : public edm::stream::EDProducer<> {
public:
  explicit FineSimFineSimTrackstersProducer(const edm::ParameterSet&);
  ~FineSimFineSimTrackstersProducer() override {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event&, const edm::EventSetup&) override;
  void computeSingleTracksterMask(std::vector<float>& result_mask,
                                  std::vector<float>& test,
                                  const Trackster& trackster,
                                  const std::vector<reco::CaloCluster>& lcs) {
    auto N = trackster.vertices().size();
    for (size_t i_lc = 0; i_lc < N; i_lc++) {
      if(lcs[trackster.vertices(i_lc)].hitsAndFractions().size() > 1){
      result_mask[trackster.vertices(i_lc)] =
          1. / trackster.vertex_multiplicity(i_lc);  // make the lcs available weighted by the vertex_multiplicity.
      test[trackster.vertices(i_lc)] = 1.;
      }
    }
  }

  // static methods for handling the global cache
  static std::unique_ptr<TrackstersCache> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(TrackstersCache*);

private:
  std::string detector_;
  bool doNose_;
  const std::string tfDnnLabel_;
  const edm::ESGetToken<TfGraphDefWrapper, TfGraphRecord> tfDnnToken_;
  const tensorflow::Session* tfSession_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTiles>> myAlgo_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>> myAlgoHFNose_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTiles>> myAlgoMIP_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>> myAlgoHFNoseMIP_;
  const edm::EDGetTokenT<std::vector<Trackster>> simtrackster_token_;
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;
  const edm::EDGetTokenT<std::vector<float>> filtered_layerclusters_mask_token_;
  const edm::EDGetTokenT<std::vector<float>> original_layerclusters_mask_token_;
  const edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> clustersTime_token_;
  edm::EDGetTokenT<TICLLayerTiles> layer_clusters_tiles_token_;
  edm::EDGetTokenT<TICLLayerTilesHFNose> layer_clusters_tiles_hfnose_token_;
  const edm::EDGetTokenT<std::vector<TICLSeedingRegion>> seeding_regions_token_;
  const std::string itername_;
  ticl::Trackster::IterationIndex iterIndex_ = ticl::Trackster::IterationIndex(0);
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geom_token_;
  hgcal::RecHitTools rhtools_;
};
DEFINE_FWK_MODULE(FineSimFineSimTrackstersProducer);

FineSimFineSimTrackstersProducer::FineSimFineSimTrackstersProducer(const edm::ParameterSet& ps)
    : detector_(ps.getParameter<std::string>("detector")),
      doNose_(detector_ == "HFNose"),
      tfDnnLabel_(ps.getParameter<std::string>("tfDnnLabel")),
      tfDnnToken_(esConsumes(edm::ESInputTag("", tfDnnLabel_))),
      tfSession_(nullptr),
      simtrackster_token_(consumes<std::vector<Trackster>>(ps.getParameter<edm::InputTag>("simTracksters"))),
      clusters_token_(consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_clusters"))),
      filtered_layerclusters_mask_token_(consumes<std::vector<float>>(ps.getParameter<edm::InputTag>("filtered_mask"))),
      original_layerclusters_mask_token_(consumes<std::vector<float>>(ps.getParameter<edm::InputTag>("original_mask"))),
      clustersTime_token_(
          consumes<edm::ValueMap<std::pair<float, float>>>(ps.getParameter<edm::InputTag>("time_layerclusters"))),
      seeding_regions_token_(
          consumes<std::vector<TICLSeedingRegion>>(ps.getParameter<edm::InputTag>("seeding_regions"))),
      itername_(ps.getParameter<std::string>("itername")),
      geom_token_(esConsumes()){
  auto plugin = ps.getParameter<std::string>("patternRecognitionBy");
  auto pluginMIP = ps.getParameter<std::string>("patternRecognitionMIPBy");
  auto pluginPSet = ps.getParameter<edm::ParameterSet>("pluginPatternRecognitionBy" + plugin);
  auto pluginPSetMIP = ps.getParameter<edm::ParameterSet>("pluginPatternRecognitionMIPBy" + pluginMIP);
  if (doNose_) {
    myAlgoHFNose_ = PatternRecognitionHFNoseFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionBy"), pluginPSet, consumesCollector());
    layer_clusters_tiles_hfnose_token_ =
        consumes<TICLLayerTilesHFNose>(ps.getParameter<edm::InputTag>("layer_clusters_hfnose_tiles"));
    myAlgoHFNoseMIP_ = PatternRecognitionHFNoseFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionMIPBy"), pluginPSetMIP, consumesCollector());
  } else {
    myAlgo_ = PatternRecognitionFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionBy"), pluginPSet, consumesCollector());
    layer_clusters_tiles_token_ = consumes<TICLLayerTiles>(ps.getParameter<edm::InputTag>("layer_clusters_tiles"));
    myAlgoMIP_ = PatternRecognitionFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionMIPBy"), pluginPSetMIP, consumesCollector());
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
  produces<std::vector<float>>("layerClustersLocalDensity");
  produces<std::vector<float>>("layerClustersRadius");
  produces<std::vector<unsigned int>>("layerClustersSize");
  produces<std::vector<unsigned int>>("layerClustersType");
//  produces<std::vector<int>>("tracksterSeeds");
}

void FineSimFineSimTrackstersProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // hgcalMultiClusters
  edm::ParameterSetDescription desc;
  desc.add<std::string>("detector", "HGCAL");
  desc.add<edm::InputTag>("simTracksters", edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<edm::InputTag>("filtered_mask", edm::InputTag("filteredLayerClusters", "iterationLabelGoesHere"));
  desc.add<edm::InputTag>("original_mask", edm::InputTag("hgcalLayerClusters", "InitialLayerClustersMask"));
  desc.add<edm::InputTag>("time_layerclusters", edm::InputTag("hgcalLayerClusters", "timeLayerCluster"));
  desc.add<edm::InputTag>("layer_clusters_tiles", edm::InputTag("ticlLayerTileProducer"));
  desc.add<edm::InputTag>("layer_clusters_hfnose_tiles", edm::InputTag("ticlLayerTileHFNose"));
  desc.add<edm::InputTag>("seeding_regions", edm::InputTag("ticlSeedingGlobal"));
  desc.add<std::string>("patternRecognitionBy", "CLUE3D");
  desc.add<std::string>("patternRecognitionMIPBy", "CA");
  desc.add<std::string>("itername", "unknown");
  desc.add<std::string>("tfDnnLabel", "tracksterSelectionTf");

  // CA Plugin
  edm::ParameterSetDescription pluginDesc;
  pluginDesc.addNode(edm::PluginDescription<PatternRecognitionFactory>("type", "CA", true));
  desc.add<edm::ParameterSetDescription>("pluginPatternRecognitionMIPByCA", pluginDesc);
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

void FineSimFineSimTrackstersProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
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
  auto clustersLocalDensity = std::make_unique<std::vector<float>>();
  auto clustersRadius = std::make_unique<std::vector<float>>();
  auto clustersSize = std::make_unique<std::vector<unsigned int>>();
  auto clustersType = std::make_unique<std::vector<unsigned int>>();
  auto output_mask = std::make_unique<std::vector<float>>();
  const std::vector<float>& original_layerclusters_mask = evt.get(original_layerclusters_mask_token_);

  tfSession_ = es.getData(tfDnnToken_).getSession();
  rhtools_.setGeometry(geom);

  // if it's regional iteration and there are seeding regions
  int tot_initial_v = 0;
  for (size_t i_st = 0; i_st < simTracksters.size(); i_st++) {
    tot_initial_v += simTracksters[i_st].vertices().size();
    int tot_number_of_lcs_in_fineSimTracksters = 0;
    int tot_number_of_lcs_in_simTracksters = 0;
    tot_number_of_lcs_in_simTracksters += simTracksters[i_st].vertices().size();
    std::vector<uint> fine_sim_trackster_index;
    std::vector<float> fine_input_cluster_mask(layerClusters.size(), 0.);
    std::vector<float> test_fine_input_cluster_mask(layerClusters.size(), 0.);
    computeSingleTracksterMask(
        fine_input_cluster_mask, test_fine_input_cluster_mask, simTracksters[i_st], layerClusters);
    std::vector<float> output_mask(fine_input_cluster_mask);
//    auto count_av = 0;
//    auto count_av_2 =
//        std::count_if(fine_input_cluster_mask.begin(), fine_input_cluster_mask.end(), [](float i) { return i > 0.; });
//    for (size_t f_i = 0; f_i < fine_input_cluster_mask.size(); ++f_i) {
//      if (fine_input_cluster_mask[f_i] > 0.) {
//        count_av += 1;
//      }
//    }
//    std::cout << "INITIAL LCs " << count_av << std::endl;
//    if (count_av != static_cast<int>(simTracksters[i_st].vertices().size()) ||
//        count_av_2 != static_cast<int>(simTracksters[i_st].vertices().size())) {
//      std::cout << "DIFFERENCE " << count_av - static_cast<int>(simTracksters[i_st].vertices().size()) << std::endl;
//      std::cout << "DIFFERENCE_2  " << count_av_2 - static_cast<int>(simTracksters[i_st].vertices().size())
//                << std::endl;
//
//      for (size_t j = 0; j < simTracksters[i_st].vertices().size(); ++j) {
//        auto lc_id = simTracksters[i_st].vertices(j);
//        if (!(fine_input_cluster_mask[lc_id] > 0.)) {
//          std::cout << " Mask " << fine_input_cluster_mask[lc_id] << " Vertex Multiplicity "
//                    << simTracksters[i_st].vertex_multiplicity(j) << std::endl;
//        }
//      }
//    }
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
                                                                                           seeding_regions,
                                                                                           tfSession_);

      typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Outputs output(
          *tmp_result, *tmp_tracksterSeeds, *clustersLocalDensity, *clustersRadius, *clustersSize, *clustersType);
      myAlgoHFNose_->makeTracksters(inputHFNose, output, seedToTrackstersAssociation);

    } else {
      const auto& layer_clusters_tiles = evt.get(layer_clusters_tiles_token_);
      const typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Inputs input(evt,
                                                                               es,
                                                                               layerClusters,
                                                                               fine_input_cluster_mask,
                                                                               layerClustersTimes,
                                                                               layer_clusters_tiles,
                                                                               seeding_regions,
                                                                               tfSession_);

      typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Outputs output(
          *tmp_result, *tmp_tracksterSeeds, *clustersLocalDensity, *clustersRadius, *clustersSize, *clustersType);
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
    for (auto& x : output_mask) {
      if (x > 0) {
        tot_still_av += 1;
      }
    }
    //std::cout << " LCs still available " << tot_still_av << std::endl;
    if (doNose_) {
      const auto& layer_clusters_hfnose_tiles_mip = evt.get(layer_clusters_tiles_hfnose_token_);
      const typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Inputs inputHFNoseMIP(
          evt,
          es,
          layerClusters,
          output_mask,
          layerClustersTimes,
          layer_clusters_hfnose_tiles_mip,
          seeding_regions,
          tfSession_);

      typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Outputs outputMIPNose(
          *tmp_resultMIP, *tmp_tracksterSeedsMIP, *clustersLocalDensity, *clustersRadius, *clustersSize, *clustersType);
      myAlgoHFNoseMIP_->makeTracksters(inputHFNoseMIP, outputMIPNose, seedToTrackstersAssociation);

    } else {
      const auto& layer_clusters_tiles_mip = evt.get(layer_clusters_tiles_token_);
      const typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Inputs inputMIP(evt,
                                                                                  es,
                                                                                  layerClusters,
                                                                                  output_mask,
                                                                                  layerClustersTimes,
                                                                                  layer_clusters_tiles_mip,
                                                                                  seeding_regions,
                                                                                  tfSession_);

      typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Outputs outputMIP(*tmp_resultMIP,
                                                                              *tmp_tracksterSeedsMIP,
                                                                              *clustersLocalDensity,
                                                                              *clustersRadius,
                                                                              *clustersSize,
                                                                              *clustersType);
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
    for(auto& tmp_t : *tmp_result){
      tmp_t.setSeed(simTracksters[i_st].seedID(), simTracksters[i_st].seedIndex());
    }
    result->insert(result->end(), tmp_result->begin(), tmp_result->end());
    resultMIP->insert(resultMIP->end(), tmp_resultMIP->begin(), tmp_resultMIP->end());
    tracksterSeedsMIP->insert(tracksterSeedsMIP->end(), tmp_tracksterSeedsMIP->begin(), tmp_tracksterSeedsMIP->end());
    tracksterSeeds->insert(tracksterSeeds->end(), tmp_tracksterSeeds->begin(), tmp_tracksterSeeds->end());

  }
// end simtracksters loop

assignPCAtoTracksters(*result,
                            layerClusters,
                            layerClustersTimes,
                            rhtools_.getPositionLayer(rhtools_.lastLayerEE(false), false).z());

// run energy regression and ID
// energyRegressionAndID(layerClusters, result);

evt.put(std::move(result), "fine");
evt.put(std::move(resultMIP), "fineMIP");
evt.put(std::move(simTracksterToFineSimTracksters), "fine");
//evt.put(std::move(tracksterSeeds), "fine");
//evt.put(std::move(tracksterSeedsMIP), "fineMIP");
evt.put(std::move(clustersLocalDensity), "layerClustersLocalDensity");
evt.put(std::move(clustersRadius), "layerClustersRadius");
evt.put(std::move(clustersSize), "layerClustersSize");
evt.put(std::move(clustersType), "layerClustersType");
//evt.put(std::move(tracksterSeeds), "tracksterSeeds");
}
