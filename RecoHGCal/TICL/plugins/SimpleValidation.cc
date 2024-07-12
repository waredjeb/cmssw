#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <variant>

#include <memory>  // unique_ptr
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/transform.h"
#include "SimCalorimetry/HGCalAssociatorProducers/interface/AssociatorTools.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "DataFormats/HGCalReco/interface/Common.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"
// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using TracksterToTracksterMap =
      ticl::AssociationMap<ticl::mapWithFractionAndScore, std::vector<ticl::Trackster>, std::vector<ticl::Trackster>>;
class SimpleValidation : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SimpleValidation(const edm::ParameterSet&);
  ~SimpleValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  float global_eff_ = 0.f;
  float global_fake_ = 0.f;
  float global_merge_ = 0.f;
  float total_recoMerged_ = 0.f;
  float total_recoCLUE3D_ = 0.f;
  float total_sim_ = 0.f;
  float total_simEff_ = 0.f;

  const edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstersMerged_token_;
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloParticles_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> simTracksters_CP_token_;
  const edm::EDGetTokenT<TracksterToTracksterMap> tsRecoToSimCP_token_;
  const edm::EDGetTokenT<TracksterToTracksterMap> tsSimToRecoCP_token_;
  const edm::EDGetTokenT<TracksterToTracksterMap> MergeRecoToSimCP_token_;
  const edm::EDGetTokenT<TracksterToTracksterMap> MergeSimToRecoCP_token_;

  TTree* output_tree_;
};

SimpleValidation::SimpleValidation(const edm::ParameterSet& iConfig)
    : tracksters_token_(
          consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstersclue3d"))),
      trackstersMerged_token_(
          consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstersMerged"))),
      caloParticles_token_(consumes<std::vector<CaloParticle>>(iConfig.getParameter<edm::InputTag>("caloParticles"))),
      simTracksters_CP_token_(
          consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("simtrackstersCP"))),
      tsRecoToSimCP_token_(consumes<TracksterToTracksterMap>(iConfig.getParameter<edm::InputTag>("recoToSimAssociatorCP"))),
      tsSimToRecoCP_token_(consumes<TracksterToTracksterMap>(iConfig.getParameter<edm::InputTag>("simToRecoAssociatorCP"))),
      MergeRecoToSimCP_token_(
          consumes<TracksterToTracksterMap>(iConfig.getParameter<edm::InputTag>("MergerecoToSimAssociatorCP"))),
      MergeSimToRecoCP_token_(
          consumes<TracksterToTracksterMap>(iConfig.getParameter<edm::InputTag>("MergesimToRecoAssociatorCP")))

{}

SimpleValidation::~SimpleValidation() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
  // if (trackLabels_[0].label().compare("pixelTracks0") == 0) {
  //   std::cerr << "pixelTracks" << "\n"
  //             << "Total Simulated "<< global_st_ << "\n"
  //             << "Total Reconstructed " << global_rt_ << "\n"
  //             << "Total Associated (recoToSim) " << global_at_ << "\n"
  //             << "Total Fakes " << global_rt_ - global_at_ << "\n"
  //             << "Total Associated (simRoReco) " << global_ast_ << "\n"
  //             << "Total Duplicated " << global_dt_ << "\n";
  // }
}

//
// member functions
//

// ------------ method called for each event  ------------
void SimpleValidation::analyze(const edm::Event& event, const edm::EventSetup& iSetup) {
  edm::Handle<std::vector<ticl::Trackster>> tracksters_handle;
  event.getByToken(tracksters_token_, tracksters_handle);
  const auto& tracksters = *tracksters_handle;
  //merged tracksters
  edm::Handle<std::vector<ticl::Trackster>> trackstersMerged_handle;
  event.getByToken(trackstersMerged_token_, trackstersMerged_handle);
  const auto& trackstersMerged = *trackstersMerged_handle;
  // simTracksters from CP
  edm::Handle<std::vector<ticl::Trackster>> simTrackstersCP_h;
  event.getByToken(simTracksters_CP_token_, simTrackstersCP_h);
  const auto& simTrackstersCP = *simTrackstersCP_h;

  // trackster reco to sim CP
  edm::Handle<TracksterToTracksterMap> tsRecoToSimCP_h;
  event.getByToken(tsRecoToSimCP_token_, tsRecoToSimCP_h);
  auto const& tsRecoSimCPMap = *tsRecoToSimCP_h;

  // sim simTrackster CP to reco trackster
  edm::Handle<TracksterToTracksterMap> tsSimToRecoCP_h;
  event.getByToken(tsSimToRecoCP_token_, tsSimToRecoCP_h);
  auto const& tsSimToRecoCPMap = *tsSimToRecoCP_h;

  edm::Handle<std::vector<CaloParticle>> caloParticles_h;
  event.getByToken(caloParticles_token_, caloParticles_h);
  auto const& caloParticles = *caloParticles_h;

  edm::Handle<TracksterToTracksterMap> mergetsRecoToSimCP_h;
  event.getByToken(MergeRecoToSimCP_token_, mergetsRecoToSimCP_h);
  auto const& MergetsRecoSimCPMap = *mergetsRecoToSimCP_h;

  // sim simTrackster CP to reco trackster
  edm::Handle<TracksterToTracksterMap> mergetsSimToRecoCP_h;
  event.getByToken(MergeSimToRecoCP_token_, mergetsSimToRecoCP_h);
  auto const& MergetsSimToRecoCPMap = *mergetsSimToRecoCP_h;


  //reco_sim_ratio = tracksters.size() / simTrackstersCP.size();

  std::vector<int> stsInTracksterC3D(tracksters.size(), 0);
  std::vector<int> stsInTracksterSignal(tracksters.size(), 0);
  std::vector<size_t> cPIndices;
  //  std::cout << "caloParticles size " << caloParticles.size() << std::endl;
  removeCPFromPU(caloParticles, cPIndices, true);
  total_sim_ += simTrackstersCP.size();
  total_simEff_ += cPIndices.size();
  total_recoCLUE3D_ += tracksters.size(); 
  total_recoMerged_ += trackstersMerged.size(); 
  //  std::cout << "cpIndices " << cPIndices.size() << std::endl;
  for (size_t iReco = 0; iReco < tracksters.size(); ++iReco) {
    // CLUE3D -> STS-CP
    const auto stsCP_vec = tsRecoSimCPMap.at(iReco);
    if (!stsCP_vec.empty()) {
      for (const auto& [sts_id, sharedEnergyAndScore] : stsCP_vec) {
        if(sharedEnergyAndScore.second <= 0.6)
          stsInTracksterC3D[iReco] += 1;
      }
    }
  }

  std::vector<int> stsInTracksterMerged(trackstersMerged.size(), 0);
  std::vector<int> stsInTracksterSignalMerged(trackstersMerged.size(), 0);
  for (size_t iReco = 0; iReco < trackstersMerged.size(); ++iReco) {
    // CLUE3D -> STS-CP
    const auto stsCP_vec = MergetsRecoSimCPMap.at(iReco);
    if (!stsCP_vec.empty()) {
      for (const auto& [sts_id, sharedEnergyAndScore] : stsCP_vec) {
        if(sharedEnergyAndScore.second <= 0.6)
          stsInTracksterMerged[iReco] += 1;
      }
    }
  }


  //Merge Rate for TrackstersMerged
 
  for (size_t iReco = 0; iReco != stsInTracksterMerged.size(); iReco++) {
    if(stsInTracksterMerged[iReco] == 0){
      global_fake_ += 1; 
    }
    if(stsInTracksterMerged[iReco] >= 2){
      global_merge_ += 1; 
    }
  }

  
  for (size_t iSim = 0; iSim != simTrackstersCP.size(); iSim++) {
    bool matchedPur = false;
    bool matchedEff = false;
    auto cpIndex = simTrackstersCP[iSim].seedIndex();
    if (std::find(cPIndices.begin(), cPIndices.end(), cpIndex) == cPIndices.end()) {
      continue;
    }
    auto recoTS_CP_vec = MergetsSimToRecoCPMap[iSim]; 
    auto maxEL = std::max_element(recoTS_CP_vec.begin(), recoTS_CP_vec.end(), [](std::pair<unsigned int, std::pair<float,float>>& p1, std::pair<unsigned int, std::pair<float,float>>& p2){
        return p1.second.first < p2.second.first;
        }
     );
    auto const maxSharedEnergy = maxEL->first;
    auto const index_maxSharedEnergy = std::distance(recoTS_CP_vec.begin(), maxEL);
    if(maxSharedEnergy / simTrackstersCP[iSim].raw_energy() >= 0.5 and matchedEff == false){
      global_eff_ += 1;
      matchedEff = true;
    }
  }

//  for (size_t iSim = 0; iSim != simTrackstersCP.size(); iSim++) {
//    bool matchedPur = false;
//    bool matchedEff = false;
//    int totMatched = 0;
//    bool merged = 0;
//    auto cpIndex = simTrackstersCP[iSim].seedIndex();
//    if(simTrackstersCP[iSim].raw_energy() >= 50.f){
//        global_reco_eff_high_energy_den += 1;
//      }
//      else{
//        global_reco_eff_low_energy_den += 1;
//      }
//      if(simTrackstersCP[iSim].barycenter().eta() >= 2.2f){
//        global_reco_eff_high_eta_den += 1;
//      }
//      else{
//        global_reco_eff_low_eta_den += 1;
//      }
//      if (std::find(cPIndices.begin(), cPIndices.end(), cpIndex) == cPIndices.end()) {
//        continue;
//      }
//    const edm::Ref<ticl::TracksterCollection> stsCPRef(simTrackstersCP_h, iSim);
//    auto const ts_iter = tsSimToRecoCPMap.find(stsCPRef);
//    if (ts_iter != tsSimToRecoCPMap.end()) {
//      const auto& tsAssociated = ts_iter->val;
//      for (auto const& ts : tsAssociated) {
//        auto ts_idx = (ts.first).get() - (edm::Ref<ticl::TracksterCollection>(tracksters_handle, 0)).get();
//        auto const& recoRef = edm::Ref<ticl::TracksterCollection>(tracksters_handle, ts_idx);
//        if (ts.second.second <= 0.2 and !matchedPur) {
//          auto const& pu_iter = tsRecoSimPUMap.find(recoRef);
//          if (pu_iter != tsRecoSimPUMap.end()) {
//            auto const& puPair = pu_iter->val[0];
//            if (puPair.second.first / recoRef->raw_energy() <= 0.1) {
//              global_ass_pur_ += 1;
//              matchedPur = true;
//            }
//          }
//        }
//        if (ts.second.first / simTrackstersCP[iSim].raw_energy() >= 0.5 and !matchedEff) {
//          // check PU contamination
//    //      if(tracksters[ts_idx].raw_energy() / simTrackstersCP[iSim].raw_energy() >= 1.0f){
//    //        global_merge_ += 1;
//    //      }
//          auto const& pu_iter = tsRecoSimPUMap.find(recoRef);
//          if (pu_iter != tsRecoSimPUMap.end()) {
//            auto const& puPair = pu_iter->val[0];
//            if (puPair.second.first / recoRef->raw_energy() <= 0.1) {
//              global_ass_eff_ += 1;
//              matchedEff = true;
//              if(simTrackstersCP[iSim].raw_energy() >= 50.f){
//                global_reco_eff_high_energy_num += 1;
//              }
//              else{
//                global_reco_eff_low_energy_num += 1;
//              }
//              if(simTrackstersCP[iSim].barycenter().eta() >= 2.2f){
//                global_reco_eff_high_eta_num += 1;
//              }
//              else{
//                global_reco_eff_low_eta_num += 1;
//              }
//            }
//          }
//        }
//      }
//    }
//  }
//  global_reco_ += tracksters.size();
//  global_sim_ += simTrackstersCP.size();
}

// ------------ method called once each job just before starting event loop  ------------
void SimpleValidation::beginJob() {
  // please remove this method if not needed
  edm::Service<TFileService> fs;
  output_tree_ = fs->make<TTree>("output", "putput params");

  output_tree_->Branch("total_recoCLUE3D", &total_recoCLUE3D_);
  output_tree_->Branch("total_recoMerged", &total_recoMerged);
  output_tree_->Branch("number_of_fake", &global_fake_);
  output_tree_->Branch("number_of_merge", &global_merge_);
  output_tree_->Branch("number_of_sim", &total_sim_);
  output_tree_->Branch("number_of_sim_eff", &total_simEff_);
  output_tree_->Branch("number_of_eff", &global_eff_);
}

// ------------ method called once each job just after ending the event loop  ------------
void SimpleValidation::endJob() {
  // please remove this method if not needed
  output_tree_->Fill();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SimpleValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackstersclue3d", edm::InputTag("ticlTrackstersCLUE3D"));
  desc.add<edm::InputTag>("trackstersMerged", edm::InputTag("ticlCandidate"));
  desc.add<edm::InputTag>("simtrackstersCP", edm::InputTag("ticlSimTracksters", "fromCPs"));
  desc.add<edm::InputTag>("caloParticles", edm::InputTag("mix", "MergedCaloTruth"));
  //  desc.add<edm::InputTag>("simtrackstersPU", edm::InputTag("ticlSimTracksters", "PU"));
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("recoToSimAssociatorCP",
                          edm::InputTag("tracksterSimTracksterFromCPsAssociationPR", "tracksterToSimTracksterMap"));
  desc.add<edm::InputTag>("simToRecoAssociatorCP",
                          edm::InputTag("tracksterSimTracksterFromCPsAssociationPR", "simTracksterToTracksterMap"));
  desc.add<edm::InputTag>(
      "MergerecoToSimAssociatorCP",
      edm::InputTag("tracksterSimTracksterFromCPsAssociationLinking", "tracksterToSimTracksterMap"));
  desc.add<edm::InputTag>(
      "MergesimToRecoAssociatorCP",
      edm::InputTag("tracksterSimTracksterFromCPsAssociationLinking", "simTracksterToTracksterMap"));
  descriptions.add("simpleValidation", desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleValidation);
