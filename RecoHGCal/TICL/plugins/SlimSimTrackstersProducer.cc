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

#include "TrackingTools/Records/interface/TfGraphRecord.h"
#include "RecoTracker/FinalTrackSelectors/interface/TfGraphDefWrapper.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

using namespace ticl;

class SlimSimTrackstersProducer : public edm::stream::EDProducer<edm::GlobalCache<TrackstersCache>> {
public:
  explicit SlimSimTrackstersProducer(const edm::ParameterSet&, const TrackstersCache*);
  ~SlimSimTrackstersProducer() override {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event&, const edm::EventSetup&) override;
  void makeTrackstersFromSim(edm::Event& evt,
                            const edm::EventSetup& es,
                            Trackster& simTrackster,
                            const std::vector<reco::CaloCluster>& layerClusters,
                            const edm::ValueMap<std::pair<float, float>>& layerClustersTimes,
                            std::vector<float>& mask_patternRecognition,
                            const std::vector<TICLSeedingRegion>& seeding_regions,
                            std::vector<Trackster>& result,
                            std::vector<Trackster>& resultMIP,
                            Trackster& tracksterNotClustered,
                            std::vector<int>& tracksterSeeds,
                            std::vector<int>& tracksterSeedsMIP,
                            bool doNose_);
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

  int count(std::vector<float>& mask, int target){
    int count = 0;
    for(auto& x : mask){
      if(x > target){
        count += 1;
      }
    }
    return count;
  }

  void updateMask(std::vector<Trackster>& tracksters, std::vector<float>& mask){
    for (auto& trackster : tracksters) {
      // Mask the used elements, accordingly
      for (auto const v : trackster.vertices()) {
        // TODO(rovere): for the moment we mask the layer cluster completely. In
        // the future, properly compute the fraction of usage.
        mask[v] = 0.;
      }
    }
  }


  // std::vector<float> computeMaskFromLCsize(std::vector<reco::CaloCluster>& lcs, int min_hits, hgcal::RecHitTools& rhtools_){
  //   std::vector<float> output_mask = 
  //   for(size_t i = 0; i < lcs.size(); i++){
  //       auto const& layerCluster = lcs[i];
  //       auto const& haf = layerCluster.hitsAndFractions();
  //       if(haf.size() > min_hits){
          
  //       }
  //   }    
  // }
  int findNearestSeed(const reco::CaloCluster& lc,
                      std::vector<int>& tracksterSeeds,
                      const std::vector<reco::CaloCluster>& lcs) {
    double min = 1e9;
    int result = -1;
    auto distance_sq = [](reco::CaloCluster lc1, reco::CaloCluster lc2) {
      auto dx = lc2.x() - lc1.x();
      auto dy = lc2.y() - lc1.y();
      auto dz = lc2.z() - lc1.z();
      auto dx2 = dx * dx;
      auto dy2 = dy * dy;
      auto dz2 = dz * dz;

      return sqrt(dx2 + dy2 + dz2);
    };
    auto r = 0;
    for (auto& i_s : tracksterSeeds) {
      auto distance = distance_sq(lc, lcs[i_s]);
      if (distance < min) {
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
  const std::string tfDnnLabel_;
  const edm::ESGetToken<TfGraphDefWrapper, TfGraphRecord> tfDnnToken_;
  const tensorflow::Session* tfSession_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTiles>> myAlgoHigh_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>> myAlgoHFNoseHigh_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTiles>> myAlgoLow_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>> myAlgoHFNoseLow_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTiles>> myAlgoMIP_;
  std::unique_ptr<PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>> myAlgoHFNoseMIP_;

  const edm::EDGetTokenT<std::vector<Trackster>> simtrackster_token_;
  const edm::EDGetTokenT<std::vector<Trackster>> simtracksterCP_token_;
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
DEFINE_FWK_MODULE(SlimSimTrackstersProducer);

std::unique_ptr<TrackstersCache> SlimSimTrackstersProducer::initializeGlobalCache(const edm::ParameterSet& params) {
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

void SlimSimTrackstersProducer::globalEndJob(TrackstersCache* cache) {
  delete cache->eidGraphDef;
  cache->eidGraphDef = nullptr;
}

SlimSimTrackstersProducer::SlimSimTrackstersProducer(const edm::ParameterSet& ps, const TrackstersCache* cache)
    : detector_(ps.getParameter<std::string>("detector")),
      doNose_(detector_ == "HFNose"),
      tfDnnLabel_(ps.getParameter<std::string>("tfDnnLabel")),
      tfDnnToken_(esConsumes(edm::ESInputTag("", tfDnnLabel_))),
      tfSession_(nullptr),
      simtrackster_token_(consumes<std::vector<Trackster>>(ps.getParameter<edm::InputTag>("simTracksters"))),
      simtracksterCP_token_(consumes<std::vector<Trackster>>(ps.getParameter<edm::InputTag>("simTrackstersCP"))),
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
  auto pluginHigh = ps.getParameter<std::string>("patternRecognitionHighBy");
  auto pluginLow = ps.getParameter<std::string>("patternRecognitionLowBy");
  auto pluginMIP = ps.getParameter<std::string>("patternRecognitionMIPBy");
  auto pluginPSetHigh = ps.getParameter<edm::ParameterSet>("pluginPatternRecognitionHighBy" + pluginHigh);
  auto pluginPSetLow = ps.getParameter<edm::ParameterSet>("pluginPatternRecognitionLowBy" + pluginLow);
  auto pluginPSetMIP = ps.getParameter<edm::ParameterSet>("pluginPatternRecognitionMIPBy" + pluginMIP);
  if (doNose_) {
    myAlgoHFNoseHigh_ = PatternRecognitionHFNoseFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionHighBy"), pluginPSetHigh, consumesCollector());
    myAlgoHFNoseLow_ = PatternRecognitionHFNoseFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionLowBy"), pluginPSetLow, consumesCollector());
    layer_clusters_tiles_hfnose_token_ =
        consumes<TICLLayerTilesHFNose>(ps.getParameter<edm::InputTag>("layer_clusters_hfnose_tiles"));
    myAlgoHFNoseMIP_ = PatternRecognitionHFNoseFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionMIPBy"), pluginPSetMIP, consumesCollector());
  } else {
    myAlgoHigh_ = PatternRecognitionFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionHighBy"), pluginPSetHigh, consumesCollector());
    myAlgoLow_ = PatternRecognitionFactory::get()->create(
        ps.getParameter<std::string>("patternRecognitionLowBy"), pluginPSetLow, consumesCollector());
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

  produces<std::vector<Trackster>>();
  produces<std::vector<int>>();
  produces<std::vector<Trackster>>("notClustered");
  produces<std::vector<Trackster>>("notClusteredCP");
  produces<std::vector<Trackster>>("fromCPs");
  produces<std::vector<int>>("fromCPs");
  produces<std::vector<Trackster>>("slimMIP");//temporary
  produces<std::vector<int>>("slimMIP");//temporary
  produces<std::vector<Trackster>>("slimMIPCP");//temporary
  produces<std::vector<int>>("slimMIPCP");//temporary
  produces<std::map<uint, std::vector<uint>>>();
  produces<std::vector<float>>();  //  Mask to be applied at the next iteration
}

void SlimSimTrackstersProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // hgcalMultiClusters
  edm::ParameterSetDescription desc;
  desc.add<std::string>("detector", "HGCAL");
  desc.add<edm::InputTag>("simTracksters", edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>("simTrackstersCP", edm::InputTag("ticlSimTracksters", "fromCPs"));
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<edm::InputTag>("filtered_mask", edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>("time_layerclusters", edm::InputTag("hgcalLayerClusters", "timeLayerCluster"));
  desc.add<edm::InputTag>("original_mask", edm::InputTag("hgcalLayerClusters", "InitialLayerClustersMask"));
  desc.add<edm::InputTag>("trkparticles", edm::InputTag("prunedTrackingParticles"));
  desc.add<edm::InputTag>("seeding_regions", edm::InputTag("ticlSeedingGlobal"));
  desc.add<edm::InputTag>("layer_clusters_tiles", edm::InputTag("ticlLayerTileProducer"));
  desc.add<edm::InputTag>("layer_clusters_hfnose_tiles", edm::InputTag("ticlLayerTileHFNose"));
  desc.add<std::string>("patternRecognitionHighBy", "CLUE3D");
  desc.add<std::string>("patternRecognitionLowBy", "CLUE3D");
  desc.add<std::string>("patternRecognitionMIPBy", "CA");
  desc.add<std::string>("eid_graph_path", "RecoHGCal/TICL/data/tf_models/energy_id_v0.pb");
  desc.add<std::string>("itername", "unknown");
  desc.add<std::string>("tfDnnLabel", "tracksterSelectionTf");

  // CA Plugin
  edm::ParameterSetDescription pluginDesc;
  pluginDesc.addNode(edm::PluginDescription<PatternRecognitionFactory>("type", "CA", true));
  desc.add<edm::ParameterSetDescription>("pluginPatternRecognitionMIPByCA", pluginDesc);
  //CA Plugin for MIP
  // edm::ParameterSetDescription pluginDescCA;
  // pluginDescCA.addNode(edm::PluginDescription<PatternRecognitionFactory>("type", "CA", true));
  // desc.add<edm::ParameterSetDescription>("pluginPatternRecognitionByCA", pluginDescCA);
  // CLUE3D Plugin
  edm::ParameterSetDescription pluginDescClue3DHigh;
  pluginDescClue3DHigh.addNode(edm::PluginDescription<PatternRecognitionFactory>("type", "CLUE3D", true));
  desc.add<edm::ParameterSetDescription>("pluginPatternRecognitionHighByCLUE3D", pluginDescClue3DHigh);

  edm::ParameterSetDescription pluginDescClue3DLow;
  pluginDescClue3DLow.addNode(edm::PluginDescription<PatternRecognitionFactory>("type", "CLUE3D", true));
  desc.add<edm::ParameterSetDescription>("pluginPatternRecognitionLowByCLUE3D", pluginDescClue3DLow);


  descriptions.add("slimSimTrackstersProducer", desc);
}

void SlimSimTrackstersProducer::makeTrackstersFromSim(edm::Event& evt,
                                                      const edm::EventSetup& es,
                                                      Trackster& simTrackster,
                                                      const std::vector<reco::CaloCluster>& layerClusters,
                                                      const edm::ValueMap<std::pair<float, float>>& layerClustersTimes,
                                                      std::vector<float>& mask_patternRecognition,
                                                      const std::vector<TICLSeedingRegion>& seeding_regions,
                                                      std::vector<Trackster>& result,
                                                      std::vector<Trackster>& resultMIP,
                                                      Trackster& tracksterFromNotClusteredLCs,
                                                      std::vector<int>& tracksterSeeds,
                                                      std::vector<int>& tracksterSeedsMIP,
                                                      bool doNose_) {
  auto tmp_result = std::make_unique<std::vector<Trackster>>();
  auto tmp_tracksterSeeds = std::make_unique<std::vector<int>>();
  auto tmp_resultMIP = std::make_unique<std::vector<Trackster>>();
  auto tmp_tracksterSeedsMIP = std::make_unique<std::vector<int>>();
  auto tmp_tracksterSeedsDoubletsMIP = std::make_unique<std::vector<std::vector<int>>>();
  auto tmp_tracksterSeedsDoublets = std::make_unique<std::vector<std::vector<int>>>();
  // Trackster tracksterFromNotClusteredLCs();
  std::vector<float> output_mask(mask_patternRecognition);
  std::unordered_map<int, std::vector<int>> seedToTrackstersAssociation;
  
  // if it's regional iteration and there are seeding regions
  if (!seeding_regions.empty() and seeding_regions[0].index != -1) {
    auto numberOfSeedingRegions = seeding_regions.size();
    for (unsigned int i = 0; i < numberOfSeedingRegions; ++i) {
      seedToTrackstersAssociation.emplace(seeding_regions[i].index, 0);
    }
  }
  tfSession_ = es.getData(tfDnnToken_).getSession();
  if (doNose_) {
    const auto& layer_clusters_hfnose_tiles = evt.get(layer_clusters_tiles_hfnose_token_);
    const typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Inputs inputHFNose(evt,
                                                                                         es,
                                                                                         layerClusters,
                                                                                         mask_patternRecognition,
                                                                                         layerClustersTimes,
                                                                                         layer_clusters_hfnose_tiles,
                                                                                         seeding_regions,
                                                                                         tfSession_);

    typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Outputs output(
        *tmp_result, *tmp_tracksterSeeds, *tmp_tracksterSeedsDoublets);
    myAlgoHFNoseHigh_->makeTracksters(inputHFNose, output, seedToTrackstersAssociation);

  } else {
    const auto& layer_clusters_tiles = evt.get(layer_clusters_tiles_token_);
    const typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Inputs input(
        evt, es, layerClusters, mask_patternRecognition, layerClustersTimes, layer_clusters_tiles, seeding_regions, tfSession_);

    typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Outputs output(
        *tmp_result, *tmp_tracksterSeeds, *tmp_tracksterSeedsDoublets);
    myAlgoHigh_->makeTracksters(input, output, seedToTrackstersAssociation);
  }
  
  updateMask(*tmp_result, output_mask);
  

  // if (doNose_) {
  //   const auto& layer_clusters_hfnose_tiles = evt.get(layer_clusters_tiles_hfnose_token_);
  //   const typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Inputs inputHFNose(evt,
  //                                                                                        es,
  //                                                                                        layerClusters,
  //                                                                                        output_mask,
  //                                                                                        layerClustersTimes,
  //                                                                                        layer_clusters_hfnose_tiles,
  //                                                                                        seeding_regions,
  //                                                                                        tfSession_);

  //   typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Outputs output(
  //       *tmp_result, *tmp_tracksterSeeds, *tmp_tracksterSeedsDoublets);
  //   myAlgoHFNoseLow_->makeTracksters(inputHFNose, output, seedToTrackstersAssociation);

  // } else {
  //   const auto& layer_clusters_tiles = evt.get(layer_clusters_tiles_token_);
  //   const typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Inputs input(
  //       evt, es, layerClusters, output_mask, layerClustersTimes, layer_clusters_tiles, seeding_regions, tfSession_);

  //   typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Outputs output(
  //       *tmp_result, *tmp_tracksterSeeds, *tmp_tracksterSeedsDoublets);
  //   myAlgoLow_->makeTracksters(input, output, seedToTrackstersAssociation);
  // }

  if (tmp_tracksterSeeds->empty()) {
    tmp_result->push_back(simTrackster);
  } 

  // std::cout << "BEFORE MASK FOR CA AVAILABLE" << count(output_mask, 0) << std::endl;
  

  // LCs recovery
  //if no seeds, just take the simtrackster.
  if (tmp_tracksterSeeds->empty()) {
    tmp_result->push_back(simTrackster);
  }
  
  updateMask(*tmp_result, output_mask);
  //RUN CA on missing LCs
  
  if (doNose_) {
    const auto& layer_clusters_hfnose_tiles_mip = evt.get(layer_clusters_tiles_hfnose_token_);
    const typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Inputs inputHFNoseMIP(
        evt, es, layerClusters, output_mask, layerClustersTimes, layer_clusters_hfnose_tiles_mip, seeding_regions, tfSession_);

    typename PatternRecognitionAlgoBaseT<TICLLayerTilesHFNose>::Outputs outputMIPNose(
        *tmp_resultMIP, *tmp_tracksterSeedsMIP, *tmp_tracksterSeedsDoubletsMIP);
    myAlgoHFNoseMIP_->makeTracksters(inputHFNoseMIP, outputMIPNose, seedToTrackstersAssociation);

  } else {
    const auto& layer_clusters_tiles_mip = evt.get(layer_clusters_tiles_token_);
    const typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Inputs inputMIP(
        evt, es, layerClusters, output_mask, layerClustersTimes, layer_clusters_tiles_mip, seeding_regions, tfSession_);

    typename PatternRecognitionAlgoBaseT<TICLLayerTiles>::Outputs outputMIP(
        *tmp_resultMIP, *tmp_tracksterSeedsMIP, *tmp_tracksterSeedsDoubletsMIP);
    myAlgoMIP_->makeTracksters(inputMIP, outputMIP, seedToTrackstersAssociation);
  }


  updateMask(*tmp_resultMIP, output_mask);

  for (size_t i = 0; i < simTrackster.vertices().size(); i++) {
    auto i_lc_st = simTrackster.vertices(i);
    auto index_nearest_seed = -1;
    if (output_mask[i_lc_st] > 0.) {
      auto lc = layerClusters[i_lc_st];
      tracksterFromNotClusteredLCs.vertices().push_back(i_lc_st);
      index_nearest_seed = findNearestSeed(lc, *tmp_tracksterSeeds, layerClusters);
      edm::LogVerbatim("SlimSimTrackstersProducer") <<  "Index nearest seed " << index_nearest_seed << " Mask " << output_mask[i_lc_st] << std::endl;
      if (index_nearest_seed >= 0) {
        (*tmp_result)[index_nearest_seed].vertices().push_back(i_lc_st);
      }
    }
  }

  result.insert(result.end(), tmp_result->begin(), tmp_result->end());
  // result.insert(result.end(), tmp_resultMIP->begin(), tmp_resultMIP->end());
  resultMIP.insert(resultMIP.end(), tmp_resultMIP->begin(), tmp_resultMIP->end());
  tracksterSeedsMIP.insert(tracksterSeedsMIP.end(), tmp_tracksterSeedsMIP->begin(), tmp_tracksterSeedsMIP->end());
  tracksterSeeds.insert(tracksterSeeds.end(), tmp_tracksterSeeds->begin(), tmp_tracksterSeeds->end());

  tmp_result->clear();
  tmp_tracksterSeeds->clear();
  
  for (auto& t : result) {
    t.setSeed(simTrackster.seedID(), simTrackster.seedIndex());
    t.setIteration(ticl::Trackster::SIM);
    t.setProbabilities(const_cast<float*>(&(simTrackster.id_probabilities()[0])));
    auto energy = 0.;
    for (size_t i_lc = 0; i_lc != t.vertices().size(); i_lc++) {
      energy += layerClusters[t.vertices(i_lc)].energy();  // computing raw energy
    }
    t.setRegressedEnergy(energy);
    t.setRawEnergy(energy);
  }
  
  // for (auto& t : trackster) {
    tracksterFromNotClusteredLCs.setSeed(simTrackster.seedID(), simTrackster.seedIndex());
    tracksterFromNotClusteredLCs.setIteration(ticl::Trackster::SIM);
    tracksterFromNotClusteredLCs.setProbabilities(const_cast<float*>(&(simTrackster.id_probabilities()[0])));
    auto energy = 0.;
    for (size_t i_lc = 0; i_lc != tracksterFromNotClusteredLCs.vertices().size(); i_lc++) {
      energy += layerClusters[tracksterFromNotClusteredLCs.vertices(i_lc)].energy();  // computing raw energy
    }
    tracksterFromNotClusteredLCs.setRegressedEnergy(energy);
    tracksterFromNotClusteredLCs.setRawEnergy(energy);
  // }
}  // end simtracksters loop

void SlimSimTrackstersProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
  std::cout << " ***** Begin SlimSimTrackstersProducer *****" << std::endl;
  auto result = std::make_unique<std::vector<Trackster>>();
  auto resultNotClustered = std::make_unique<std::vector<Trackster>>();
  auto resultMIP = std::make_unique<std::vector<Trackster>>();
  auto resultNotClusteredCP = std::make_unique<std::vector<Trackster>>();
  auto tracksterSeedsMIP = std::make_unique<std::vector<int>>();
  auto tracksterSeeds = std::make_unique<std::vector<int>>();
  auto tracksterSeedsDoublets = std::make_unique<std::vector<std::vector<int>>>();

  auto resultCP = std::make_unique<std::vector<Trackster>>();
  auto resultMIPCP = std::make_unique<std::vector<Trackster>>();
  auto tracksterSeedsMIPCP = std::make_unique<std::vector<int>>();
  auto tracksterSeedsCP = std::make_unique<std::vector<int>>();
  auto tracksterSeedsDoubletsCP = std::make_unique<std::vector<std::vector<int>>>();

  Trackster tracksterNotClusteredCP;
  Trackster tracksterNotClustered;
  auto layer_clusters_tiles = std::make_unique<TICLLayerTiles>();
  auto layer_clusters_hfnose_tiles = std::make_unique<TICLLayerTilesHFNose>();
  auto simTracksterToFineSimTracksters = std::make_unique<std::map<uint, std::vector<uint>>>();
  
  const std::vector<Trackster>& simTracksters = evt.get(simtrackster_token_);
  const std::vector<Trackster>& simTrackstersCP = evt.get(simtracksterCP_token_);
  const auto& layerClusters = evt.get(clusters_token_);
  const auto& layerClustersTimes = evt.get(clustersTime_token_);
  const auto& seeding_regions = evt.get(seeding_regions_token_);
  const auto& geom = es.getData(geom_token_);
  rhtools_.setGeometry(geom);

  // loop over simtracksters from caloparticle
  for (size_t i_st = 0; i_st < simTrackstersCP.size(); i_st++) {
    int tot_number_of_lcs_in_fineSimTracksters = 0;
    int tot_number_of_lcs_in_simTracksters = 0;
    auto simTracksterCP = simTrackstersCP[i_st];
    tot_number_of_lcs_in_simTracksters += simTracksterCP.vertices().size();
    std::vector<float> fine_input_cluster_mask(layerClusters.size(), 0.);
    std::vector<float> test_fine_input_cluster_mask(layerClusters.size(), 0.);
    computeSingleTracksterMask(fine_input_cluster_mask, test_fine_input_cluster_mask, simTracksterCP, layerClusters);
    makeTrackstersFromSim(evt,
                          es,
                          simTracksterCP,
                          layerClusters,
                          layerClustersTimes,
                          fine_input_cluster_mask,
                          seeding_regions,
                          *resultCP,
                          *resultMIPCP,
                          tracksterNotClusteredCP,
                          *tracksterSeedsCP,
                          *tracksterSeedsMIPCP,
                          doNose_);
  }

  for (size_t i_st = 0; i_st < simTracksters.size(); i_st++) {
    std::cout << "FROM SIM CLUSTER " << i_st << std::endl;
    std::vector<uint> fine_sim_trackster_index;
    int tot_number_of_lcs_in_fineSimTracksters = 0;
    int tot_number_of_lcs_in_simTracksters = 0;
    auto simTrackster = simTracksters[i_st];
    
    tot_number_of_lcs_in_simTracksters += simTrackster.vertices().size();
    std::vector<float> fine_input_cluster_mask(layerClusters.size(), 0.);
    std::vector<float> test_fine_input_cluster_mask(layerClusters.size(), 0.);
    computeSingleTracksterMask(fine_input_cluster_mask, test_fine_input_cluster_mask, simTrackster, layerClusters);
    makeTrackstersFromSim(evt,
                          es,
                          simTrackster,
                          layerClusters,
                          layerClustersTimes,
                          fine_input_cluster_mask,
                          seeding_regions,
                          *result,
                          *resultMIP,
                          tracksterNotClustered,
                          *tracksterSeeds,
                          *tracksterSeedsMIP,
                          doNose_);
  for(size_t ind = 0; ind < result->size(); ind++){
    fine_sim_trackster_index.push_back(ind);
  }
  (*simTracksterToFineSimTracksters)[i_st] = fine_sim_trackster_index;
  }

  std::cout << " Number of SlimSimTrackster " << result->size() << " FROM MIP " << resultMIP->size() << std::endl;
  std::cout << " Number of SlimSimTrackster from CP " << resultCP->size() << " FROM MIP " << resultMIPCP->size() << std::endl;
  std::cout << " Number of TracksterSeeds " << tracksterSeeds->size() << " FROM MIP " <<  tracksterSeedsMIP->size() << std::endl;
  std::cout << " Number of TracksterSeeds from CP " << tracksterSeedsCP->size() << " FROM MIP " <<  tracksterSeedsMIPCP->size() << std::endl;
  resultNotClustered->push_back(tracksterNotClustered);
  resultNotClusteredCP->push_back(tracksterNotClusteredCP);
  evt.put(std::move(result));
  evt.put(std::move(resultNotClustered),"notClustered");
  evt.put(std::move(resultNotClusteredCP),"notClusteredCP");
  evt.put(std::move(resultMIP), "slimMIP"); //temporary
  evt.put(std::move(resultCP), "fromCPs");
  evt.put(std::move(resultMIPCP), "slimMIPCP"); //temporary
  evt.put(std::move(simTracksterToFineSimTracksters));
  evt.put(std::move(tracksterSeeds));
  // evt.put(std::move(tracksterSeedsMIP), "slimMIP"); //temporary
  // evt.put(std::move(tracksterSeedsCP), "fromCPs");
  // evt.put(std::move(tracksterSeedsMIPCP), "slimMIPCP"); //temporary
  std::cout << " ***** End SlimSimTrackstersProducer *****" << std::endl;
}
