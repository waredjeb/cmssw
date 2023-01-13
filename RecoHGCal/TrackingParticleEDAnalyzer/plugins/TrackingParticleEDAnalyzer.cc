
// -*- C++ -*-
//
// Package:    RecoHGCal/TrackingParticleEDAnalyzer
// Class:      TrackingParticleEDAnalyzer
//
/**\class TrackingParticleEDAnalyzer TrackingParticleEDAnalyzer.cc RecoHGCal/TrackingParticleEDAnalyzer/plugins/TrackingParticleEDAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Wahid Redjeb
//         Created:  Mon, 07 Nov 2022 16:46:30 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/Associations/interface/TrackAssociation.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/Track/interface/UniqueSimTrackId.h"
using namespace ticl;
using namespace reco;
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class TrackingParticleEDAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TrackingParticleEDAnalyzer(const edm::ParameterSet&);
  ~TrackingParticleEDAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<reco::Track>> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<TrackingParticle>>
      trackingParticleToken_;                              //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<SimTrack>> simTrackToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<SimVertex>> simVertexToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<ticl::TICLCandidate>>
      ticlCandidatesToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<Trackster>>
      simTrackstersToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<CaloParticle>>
      caloParticlesToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<SimCluster>>
      simClustersToken_;  //used to select what tracks to read from configuration file
                          //  std::vector<edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator>> associatorTokens;
  edm::EDGetTokenT<reco::SimToRecoCollection> associatormapStRsToken_;
  edm::EDGetTokenT<reco::RecoToSimCollection> associatormapRtSsToken_;
  edm::EDGetTokenT<SimTrackToTPMap> associationSimTrackToTPToken_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrackingParticleEDAnalyzer::TrackingParticleEDAnalyzer(const edm::ParameterSet& iConfig)
    : tracksToken_(consumes<edm::View<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      trackingParticleToken_(consumes<std::vector<TrackingParticle>>(iConfig.getParameter<edm::InputTag>("tps"))),
      simTrackToken_(consumes<std::vector<SimTrack>>(iConfig.getParameter<edm::InputTag>("simTracks"))),
      simVertexToken_(consumes<std::vector<SimVertex>>(iConfig.getParameter<edm::InputTag>("simVertices"))),
      ticlCandidatesToken_(consumes<std::vector<ticl::TICLCandidate>>(iConfig.getParameter<edm::InputTag>("ticlCandidates"))),
      simTrackstersToken_(consumes<std::vector<Trackster>>(iConfig.getParameter<edm::InputTag>("simTracksters"))),
      caloParticlesToken_(consumes<std::vector<CaloParticle>>(iConfig.getParameter<edm::InputTag>("caloParticles"))),
      simClustersToken_(consumes(iConfig.getParameter<edm::InputTag>("simClusters"))),
      associatormapStRsToken_(consumes(iConfig.getParameter<edm::InputTag>("tpToTrack"))),
      associatormapRtSsToken_(consumes(iConfig.getParameter<edm::InputTag>("trackToTp"))),
      associationSimTrackToTPToken_(consumes(iConfig.getParameter<edm::InputTag>("simTrackToTPMap"))) {
  //now do what ever initialization is needed
}

TrackingParticleEDAnalyzer::~TrackingParticleEDAnalyzer() {
}

//
// member functions
//

// ------------ method called for each event  ------------
void TrackingParticleEDAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  std::cout << "Event " << iEvent.id() << std::endl;
  edm::Handle<std::vector<TrackingParticle>> trackingParticles_h;
  iEvent.getByToken(trackingParticleToken_, trackingParticles_h);
  const auto& trackingParticles = *trackingParticles_h;

  edm::Handle<edm::View<reco::Track>> tracks_h;
  iEvent.getByToken(tracksToken_, tracks_h);
  const auto& tracks = *tracks_h;

  edm::Handle<std::vector<SimTrack>> simTracks_h;
  iEvent.getByToken(simTrackToken_, simTracks_h);
  const auto& simTracks = *simTracks_h;

  edm::Handle<std::vector<SimVertex>> simVertices_h;
  iEvent.getByToken(simVertexToken_, simVertices_h);
  const auto& simVertices = *simVertices_h;

  edm::Handle<std::vector<ticl::Trackster>> simTracksters_h;
  iEvent.getByToken(simTrackstersToken_, simTracksters_h);
  const auto& simTracksters = *simTracksters_h;

  edm::Handle<std::vector<CaloParticle>> caloParticle_h;
  iEvent.getByToken(caloParticlesToken_, caloParticle_h);
  const auto& caloParticles = *caloParticle_h;

  edm::Handle<std::vector<SimCluster>> simClusters_h;
  iEvent.getByToken(simClustersToken_, simClusters_h);
  const auto& simClusters = *simClusters_h;

  edm::Handle<RecoToSimCollection> associatorRecoToSim_h;
  iEvent.getByToken(associatormapRtSsToken_, associatorRecoToSim_h);
  const auto& associatorRecoToSim = *associatorRecoToSim_h;

  edm::Handle<SimToRecoCollection> associatorSimToReco_h;
  iEvent.getByToken(associatormapStRsToken_, associatorSimToReco_h);
  const auto& associatorSimToReco = *associatorSimToReco_h;

  edm::Handle<SimTrackToTPMap> associatorSimTrackToTP_h;
  iEvent.getByToken(associationSimTrackToTPToken_, associatorSimTrackToTP_h);
  const auto& associatorSimTrackToTP = *associatorSimTrackToTP_h;

  edm::Handle<ticl::TICLCandidate> ticlCandidates_h;
  iEvent.getByToken(ticlCandidatesToken_, ticlCandidates_h);
  const auto& ticlCandidates = *ticlCandidates_h;

  std::cout << " NUMBER OF RECO TRACKS " << tracks.size() << std::endl;
  std::cout << " NUMBER OF RECO TICL Candidates " << ticlCandidates.size() << std::endl;
  std::cout << " NUMBER OF TPS " << trackingParticles.size() << std::endl;
  std::cout << " NUMBER OF SIMTRACKS " << simTracks.size() << std::endl;
  std::cout << " NUMBER OF SIMTRACKSTERS " << simTracksters.size() << std::endl;
  std::cout << " NUMBER OF CALOPARTICLES " << caloParticles.size() << std::endl;
  std::cout << "##### ASSOCIATIONS MAPS #####" << std::endl;

  for(size_t i_st = 0; i_st < simTracksters.size(); i_st++){
    auto const& simTrackster = simTracksters[i_st];
  }
//  int recoTrackMatch = 0;
//  for (size_t i = 0; i < tracks.size(); i++) {
//    auto trackRef = tracks.refAt(i);
//    const auto recoToSim_iter = associatorRecoToSim_h.product()->find(trackRef);
//
//    if (recoToSim_iter != associatorRecoToSim.end()) {
//      const auto& tpAssociated = recoToSim_iter->val;
//      //std::cout << "TP Associated " << tpAssociated.size() << " TP Available " << trackingParticles.size() << std::endl;
//      int ntpAss = 0;
//      for (auto& tpAss : tpAssociated) {
//        ntpAss++;
//        auto tp_id = (tpAss.first).get() - (edm::Ref<std::vector<TrackingParticle>>(trackingParticles_h, 0)).get();
//        auto quality = tpAss.second;
//        recoTrackMatch++;
//      }
//    } else {
//    }
//  }
//  std::cout << "Number of Matched Reco Track " << recoTrackMatch << std::endl;
//  int ats = 0;
//  for (size_t i_tp = 0; i_tp < trackingParticles.size(); i_tp++) {
//    TrackingParticleRef tpr(trackingParticles_h, i_tp);
//    auto tp = *tpr;
//
//    if (associatorSimToReco.find(tpr) != associatorSimToReco.end()) {
//      auto const& rt = associatorSimToReco[tpr];
//      if (!rt.empty()) {
//        ats++;  //This counter counts the number of simTracks that have a recoTrack associated
//      }
//    }
//  }
//  std::cout << "TP associated to a reco track " << ats << std::endl;
//    // SIM TO RECO
//    //loop over SimTracksters
//    int simTracksterMatched = 0;
//    int tracksterPassed = 0;
//    for (size_t i = 0; i < simTracksters.size(); i++) {
//      std::cout << "SimTrackster " << i << std::endl;
//      auto const& simTrackster = simTracksters[i];
//
//    for (auto const &p : simTrackster.id_probabilities()) {
//      std::cout << std::fixed << p << " ";
//    }
//    std::cout << "\n";
//      //get Simclusters from SeedIndex of SimTrackster
//      auto const& sc = simClusters[simTrackster.seedIndex()];
//      auto const& scG4Track = sc.g4Tracks()[0];
//      UniqueSimTrackId simTkIds(scG4Track.trackId(), scG4Track.eventId());
//      auto ipos = associatorSimTrackToTP.mapping.find(simTkIds);
//      if (ipos != associatorSimTrackToTP.mapping.end()) {
//        
//        auto tpIdx =  (ipos->second).get() - (edm::Ref<std::vector<TrackingParticle>>(trackingParticles_h, 0)).get();
//        std::cout << "TPIDX Found " << tpIdx << std::endl;
//        TrackingParticleRef tpRef(trackingParticles_h, tpIdx);
//        const auto simToReco_iter = associatorSimToReco.find(tpRef);
//        if(simToReco_iter != associatorSimToReco.end()){
//          const auto& tracksAssociated = simToReco_iter->val;
//          for(auto const& trackAss : tracksAssociated){
//            auto track_id = (trackAss.first).get() - (edm::Ref<edm::View<reco::Track>>(tracks_h,0)).get();
//            std::cout << " Matched To RecoTrack " << track_id << std::endl;
//          }
//        }
//        else{
//          std::cout <<" Not Matched " << std::endl;
//        }
//      }
//    }
////      auto const& simTrackster = simTracksters[i];
////      //skip neutral particles
////      if(simTrackster.id_probability(ticl::Trackster::ParticleType::photon) > 0.5 || simTrackster.id_probability(ticl::Trackster::ParticleType::neutral_pion) > 0.5 || simTrackster.id_probability(ticl::Trackster::ParticleType::neutral_hadron)){
////        continue;
////      }
////      tracksterPassed++;
////      //get SimClusters from SeedIndex of SimTrackster
////      auto const& sc = simClusters[simTrackster.seedIndex()];
////      //get associated g4Tracks associated to SimCluster
////      auto const& scG4Track = sc.g4Tracks()[0];
////      // look for TrackingParticle that contains the simTrack
////      int tpIndex = -99;
////      for (int j = 0; j < static_cast<int>(trackingParticles.size()); j++) {
////        for (auto const& simTrackTP : trackingParticles[j].g4Tracks()) {
////            if(simTrackTP.trackId() == scG4Track.trackId()) {
////              tpIndex = j;
////            }
////  //        if (scG4Track.vertIndex() == -1) {
////  //          continue;
////  //        }
////  //        //use simVertex to see if the SimTrack belongs to the TP
////  //        auto& sVsc = simVertices[scG4Track.vertIndex()].position();
////  //        if (simTrackTP.vertIndex() == -1) {
////  //          continue;
////  //        }
////  //        auto& sVtp = simVertices[simTrackTP.vertIndex()].position();
////  //        if (sVsc == sVtp) {
////  //          std::cout << "\tTRACK ID SC " << scG4Track.trackId() << " TRACK ID TP " << simTrackTP.trackId() << std::endl;
////  //          tpIndex = j;
////  //        }
////  //        else{
////  //          std::cout << "NO TRACK ID SC " << scG4Track.trackId() << " TRACK ID TP " << simTrackTP.trackId() << std::endl;
////  //        }
////        }
////      }
////      if (tpIndex != -99) {
////        auto const& tp = trackingParticles[tpIndex];
////        TrackingParticleRef tpRef(trackingParticles_h, tpIndex);
////        const auto simToReco_iter = associatorSimToReco.find(tpRef);
////        if(simToReco_iter != associatorSimToReco.end()){
////          const auto& tracksAssociated = simToReco_iter->val;
////  //        std::cout << tracksAssociated << std::endl; std::vector<std::pair<std::RefToBase<reco::Track>, double>>
////          if(!tracksAssociated.empty()){
////            simTracksterMatched += 1;
////          }
////          for(auto const& trackAss : tracksAssociated){
////            auto track_id = (trackAss.first).get() - (edm::Ref<edm::View<reco::Track>>(tracks_h, 0)).get();
////            break;
////          }
////        }else{
////        }
////      }
////    }
////    std::cout << "Sim Tracksters matched to reco track " << simTracksterMatched << std::endl;
////    std::cout << "Tracksters passed " << tracksterPassed << std::endl;
//}
//
// ------------ method called once each job just before starting event loop  ------------
void TrackingParticleEDAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void TrackingParticleEDAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TrackingParticleEDAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("tps", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("simTracks", edm::InputTag("g4SimHits"));
  desc.add<edm::InputTag>("simVertices", edm::InputTag("g4SimHits"));
  desc.add<edm::InputTag>("ticlCandidates", edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("simTracksters", edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>("caloParticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("simClusters", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("tpToTrack", edm::InputTag("trackingParticleRecoTrackAsssociation"));
  desc.add<edm::InputTag>("trackToTp", edm::InputTag("trackingParticleRecoTrackAsssociation"));
  desc.add<edm::InputTag>("simTrackToTPMap", edm::InputTag("simHitTPAssocProducer","simTrackToTP"));

  descriptions.add("trackingParticleAnalyzer", desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackingParticleEDAnalyzer);
