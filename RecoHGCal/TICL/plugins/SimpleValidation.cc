#include "TTree.h"
#include "TFile.h"
#include "TEfficiency.h"

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

#include "SimDataFormats/Associations/interface/TracksterToSimTracksterHitLCAssociator.h"
#include "RecoHGCal/TICL/interface/commons.h"
// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

class SimpleValidation : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SimpleValidation(const edm::ParameterSet&);
  ~SimpleValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  int global_reco_ = 0;
  int global_ass_pur_ = 0;
  int global_ass_eff_ = 0;
  int global_merge_ = 0;
  int global_sim_ = 0;
  int pu_reco_ = 0;
  std::vector<float> num_eff_eta;
  std::vector<float> num_eff_phi;
  std::vector<float> num_eff_energy;
  std::vector<float> num_eff_pt;

  std::vector<float> num_pur_eta;
  std::vector<float> num_pur_phi;
  std::vector<float> num_pur_energy;
  std::vector<float> num_pur_pt;

  std::vector<float> num_fake_eta;
  std::vector<float> num_fake_phi;
  std::vector<float> num_fake_energy;
  std::vector<float> num_fake_pt;

  std::vector<float> num_dup_eta;
  std::vector<float> num_dup_phi;
  std::vector<float> num_dup_energy;
  std::vector<float> num_dup_pt;

  std::vector<float> num_merge_eta;
  std::vector<float> num_merge_phi;
  std::vector<float> num_merge_energy;
  std::vector<float> num_merge_pt;

  std::vector<float> num_pu_eta;
  std::vector<float> num_pu_phi;
  std::vector<float> num_pu_energy;
  std::vector<float> num_pu_pt;

  std::vector<float> den_sim_eta;
  std::vector<float> den_sim_phi;
  std::vector<float> den_sim_energy;
  std::vector<float> den_sim_pt;

  std::vector<float> den_reco_eta;
  std::vector<float> den_reco_phi;
  std::vector<float> den_reco_energy;
  std::vector<float> den_reco_pt;

  std::vector<float> den_fake_reco_eta;
  std::vector<float> den_fake_reco_phi;
  std::vector<float> den_fake_reco_energy;
  std::vector<float> den_fake_reco_pt;


  const edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> simTracksters_CP_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> simTracksters_PU_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimCP_token_;
  const edm::EDGetTokenT<hgcal::SimToRecoCollectionSimTracksters> tsSimToRecoCP_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimPU_token_;

  TTree* output_tree_;
};

SimpleValidation::SimpleValidation(const edm::ParameterSet& iConfig)
    : tracksters_token_(
          consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstersMerged"))),
      simTracksters_CP_token_(
          consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("simtrackstersCP"))),
      simTracksters_PU_token_(
          consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("simtrackstersPU"))),
      tsRecoToSimCP_token_(consumes<hgcal::RecoToSimCollectionSimTracksters>(
          iConfig.getParameter<edm::InputTag>("recoToSimAssociatorCP"))),
      tsSimToRecoCP_token_(consumes<hgcal::SimToRecoCollectionSimTracksters>(
          iConfig.getParameter<edm::InputTag>("simToRecoAssociatorCP"))),
      tsRecoToSimPU_token_(consumes<hgcal::RecoToSimCollectionSimTracksters>(
          iConfig.getParameter<edm::InputTag>("recoToSimAssociatorPU")))

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
  // simTracksters from CP
  edm::Handle<std::vector<ticl::Trackster>> simTrackstersCP_h;
  event.getByToken(simTracksters_CP_token_, simTrackstersCP_h);
  const auto& simTrackstersCP = *simTrackstersCP_h;

  edm::Handle<std::vector<ticl::Trackster>> simTrackstersPU_h;
  event.getByToken(simTracksters_PU_token_, simTrackstersPU_h);
  const auto& simTrackstersPU = *simTrackstersPU_h;

  // trackster reco to sim CP
  edm::Handle<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimCP_h;
  event.getByToken(tsRecoToSimCP_token_, tsRecoToSimCP_h);
  auto const& tsRecoSimCPMap = *tsRecoToSimCP_h;

  // sim simTrackster CP to reco trackster
  edm::Handle<hgcal::SimToRecoCollectionSimTracksters> tsSimToRecoCP_h;
  event.getByToken(tsSimToRecoCP_token_, tsSimToRecoCP_h);
  auto const& tsSimToRecoCPMap = *tsSimToRecoCP_h;

  edm::Handle<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimPU_h;
  event.getByToken(tsRecoToSimPU_token_, tsRecoToSimPU_h);
  auto const& tsRecoSimPUMap = *tsRecoToSimPU_h;

  std::vector<int> stsInTrackster(tracksters.size(), 0);

  for (size_t iReco = 0; iReco != tracksters.size(); iReco++) {
    den_reco_pt.push_back(tracksters[iReco].raw_pt());
    den_reco_eta.push_back(tracksters[iReco].barycenter().eta());
    den_reco_phi.push_back(tracksters[iReco].barycenter().phi());
    den_reco_energy.push_back(tracksters[iReco].raw_energy());
    const edm::Ref<ticl::TracksterCollection> tsRef(tracksters_handle, iReco);
    auto const sts_iter = tsRecoSimCPMap.find(tsRef);
    if (sts_iter != tsRecoSimCPMap.end()) {
      const auto& stsAssociated = sts_iter->val;
      for (auto const& sts : stsAssociated) {
        if (sts.second.first > 0.01) {
          den_fake_reco_pt.push_back(tracksters[iReco].raw_pt());
          den_fake_reco_eta.push_back(tracksters[iReco].barycenter().eta());
          den_fake_reco_phi.push_back(tracksters[iReco].barycenter().phi());
          den_fake_reco_energy.push_back(tracksters[iReco].raw_energy());
          //				auto sts_idx = (sts.first).get() - (edm::Ref<ticl::TracksterCollection>(simTrackstersCP_h, 0)).get();
          if (sts.second.second <= 0.2) {
            stsInTrackster[iReco] += 1;
          }
        }
      }
    }
  }

  for (size_t iReco = 0; iReco != stsInTrackster.size(); iReco++) {
    if (stsInTrackster[iReco] == 1) {
      num_fake_pt.push_back(tracksters[iReco].raw_pt());
      num_fake_eta.push_back(tracksters[iReco].barycenter().eta());
      num_fake_phi.push_back(tracksters[iReco].barycenter().phi());
      num_fake_energy.push_back(tracksters[iReco].raw_energy());
    }
    else if(stsInTrackster[iReco] > 1) {
      num_dup_pt.push_back(tracksters[iReco].raw_pt());
      num_dup_eta.push_back(tracksters[iReco].barycenter().eta());
      num_dup_phi.push_back(tracksters[iReco].barycenter().phi());
      num_dup_energy.push_back(tracksters[iReco].raw_energy());
    }
  }

  for (size_t iSim = 0; iSim != simTrackstersCP.size(); iSim++) {
    bool matchedPur = false;
    bool matchedEff = false;
    den_sim_pt.push_back(simTrackstersCP[iSim].raw_pt());
    den_sim_eta.push_back(simTrackstersCP[iSim].barycenter().eta());
    den_sim_phi.push_back(simTrackstersCP[iSim].barycenter().phi());
    den_sim_energy.push_back(simTrackstersCP[iSim].raw_energy());
    int totMatched = 0;
    bool merged = 0;
    const edm::Ref<ticl::TracksterCollection> stsCPRef(simTrackstersCP_h, iSim);
    auto const ts_iter = tsSimToRecoCPMap.find(stsCPRef);
    if (ts_iter != tsSimToRecoCPMap.end()) {
      const auto& tsAssociated = ts_iter->val;
      for (auto const& ts : tsAssociated) {
        auto ts_idx = (ts.first).get() - (edm::Ref<ticl::TracksterCollection>(tracksters_handle, 0)).get();
        auto const& recoRef = edm::Ref<ticl::TracksterCollection>(tracksters_handle, ts_idx);
        if (ts.second.second <= 0.2 and !matchedPur) {
          global_ass_pur_ += 1;
          //purity
          num_pur_pt.push_back(simTrackstersCP[iSim].raw_pt());
          num_pur_eta.push_back(simTrackstersCP[iSim].barycenter().eta());
          num_pur_phi.push_back(simTrackstersCP[iSim].barycenter().phi());
          num_pur_energy.push_back(simTrackstersCP[iSim].raw_energy());
          matchedPur = true;
        }
        if (ts.second.first / simTrackstersCP[iSim].raw_energy() >= 0.5 and !matchedEff) {
          //energy eff
          num_eff_pt.push_back(simTrackstersCP[iSim].raw_pt());
          num_eff_eta.push_back(simTrackstersCP[iSim].barycenter().eta());
          num_eff_phi.push_back(simTrackstersCP[iSim].barycenter().phi());
          num_eff_energy.push_back(simTrackstersCP[iSim].raw_energy());
          matchedEff = true;
          auto const& puRef = edm::Ref<ticl::TracksterCollection>(simTrackstersPU_h, 0);
          auto const& rtsPU_iter = tsRecoSimPUMap.find(puRef);
          if (rtsPU_iter != tsRecoSimPUMap.end()) {
            int iLoop = 0;
            for (auto const& rtsPU : rtsPU_iter->val) {
              iLoop += 1;
              auto sharedEnergyPU = rtsPU.second.first;
              if ((sharedEnergyPU / tracksters[ts_idx].raw_energy()) >= 0.5f) {
                //pu contamination
                num_pu_pt.push_back(simTrackstersCP[iSim].raw_pt());
                num_pu_eta.push_back(simTrackstersCP[iSim].barycenter().eta());
                num_pu_phi.push_back(simTrackstersCP[iSim].barycenter().phi());
                num_pu_energy.push_back(simTrackstersCP[iSim].raw_energy());
                pu_reco_ += 1;
              }
            }
          }
        }
      }
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void SimpleValidation::beginJob() {
  // please remove this method if not needed
  edm::Service<TFileService> fs;
  output_tree_ = fs->make<TTree>("output", "putput params");

  output_tree_->Branch("reco", &global_reco_);
  output_tree_->Branch("recoAssPur", &global_ass_pur_);
  output_tree_->Branch("reco_ass_eff_", &global_ass_eff_);
  output_tree_->Branch("merge", &global_merge_);
  output_tree_->Branch("pu", &pu_reco_);
  output_tree_->Branch("sim", &global_sim_);
  output_tree_->Branch("num_eff_eta", &num_eff_eta);
  output_tree_->Branch("num_eff_phi", &num_eff_phi);
  output_tree_->Branch("num_eff_energy", &num_eff_energy);
  output_tree_->Branch("num_eff_pt", &num_eff_pt);

  output_tree_->Branch("num_pur_eta", &num_pur_eta);
  output_tree_->Branch("num_pur_phi", &num_pur_phi);
  output_tree_->Branch("num_pur_energy", &num_pur_energy);
  output_tree_->Branch("num_pur_pt", &num_pur_pt);

  output_tree_->Branch("num_fake_eta", &num_fake_eta);
  output_tree_->Branch("num_fake_phi", &num_fake_phi);
  output_tree_->Branch("num_fake_energy", &num_fake_energy);
  output_tree_->Branch("num_fake_pt", &num_fake_pt);

  output_tree_->Branch("num_dup_eta", &num_dup_eta);
  output_tree_->Branch("num_dup_phi", &num_dup_phi);
  output_tree_->Branch("num_dup_energy", &num_dup_energy);
  output_tree_->Branch("num_dup_pt", &num_dup_pt);

  output_tree_->Branch("num_pu_eta", &num_pu_eta);
  output_tree_->Branch("num_pu_phi", &num_pu_phi);
  output_tree_->Branch("num_pu_energy", &num_pu_energy);
  output_tree_->Branch("num_pu_pt", &num_pu_pt);

  output_tree_->Branch("den_sim_eta", &den_sim_eta);
  output_tree_->Branch("den_sim_phi", &den_sim_phi);
  output_tree_->Branch("den_sim_energy", &den_sim_energy);
  output_tree_->Branch("den_sim_pt", &den_sim_pt);

  output_tree_->Branch("den_reco_eta", &den_reco_eta);
  output_tree_->Branch("den_reco_phi", &den_reco_phi);
  output_tree_->Branch("den_reco_energy", &den_reco_energy);
  output_tree_->Branch("den_reco_pt", &den_reco_pt);

  output_tree_->Branch("den_fake_reco_eta", &den_fake_reco_eta);
  output_tree_->Branch("den_fake_reco_phi", &den_fake_reco_phi);
  output_tree_->Branch("den_fake_reco_energy", &den_fake_reco_energy);
  output_tree_->Branch("den_fake_reco_pt", &den_fake_reco_pt);

  int binsEta = 10;
  int binsEnergy = 20;
/*
  num_eff_eta = fs->make<TH1F>("num_eff_eta", "num_eff_eta", binsEta, -3.0, 3.0);
  num_eff_phi = fs->make<TH1F>("num_eff_phi", "num_eff_phi", binsEta, -3.14, 3.14);
  num_eff_energy = fs->make<TH1F>("num_eff_energy", "num_eff_energy", binsEnergy, 0.f, 300.f);
  num_eff_pt = fs->make<TH1F>("num_eff_pt", "num_eff_pt", binsEnergy, 0.f, 200.f);

  num_pur_eta = fs->make<TH1F>("num_pur_eta", "num_pur_eta", binsEta, -3.0, 3.0);
  num_pur_phi = fs->make<TH1F>("num_pur_phi", "num_pur_phi", binsEta, -3.14, 3.14);
  num_pur_energy = fs->make<TH1F>("num_pur_energy", "num_pur_energy", binsEnergy, 0.f, 300.f);
  num_pur_pt = fs->make<TH1F>("num_pur_pt", "num_pur_pt", binsEnergy, 0.f, 200.f);

  num_fake_eta = fs->make<TH1F>("num_fake_eta", "num_fake_eta", binsEta, -3.0, 3.0);
  num_fake_phi = fs->make<TH1F>("num_fake_phi", "num_fake_phi", binsEta, -3.14, 3.14);
  num_fake_energy = fs->make<TH1F>("num_fake_energy", "num_fake_energy", binsEnergy, 0.f, 300.f);
  num_fake_pt = fs->make<TH1F>("num_fake_pt", "num_fake_pt", binsEnergy, 0.f, 200.f);

  num_dup_eta = fs->make<TH1F>("num_dup_eta", "num_dup_eta", binsEta, -3.0, 3.0);
  num_dup_phi = fs->make<TH1F>("num_dup_phi", "num_dup_phi", binsEta, -3.14, 3.14);
  num_dup_energy = fs->make<TH1F>("num_dup_energy", "num_dup_energy", binsEnergy, 0.f, 300.f);
  num_dup_pt = fs->make<TH1F>("num_dup_pt", "num_dup_pt", binsEnergy, 0.f, 200.f);

  num_pu_eta = fs->make<TH1F>("num_pu_eta", "num_pu_eta", binsEta, -3.0, 3.0);
  num_pu_phi = fs->make<TH1F>("num_pu_phi", "num_pu_phi", binsEta, -3.14, 3.14);
  num_pu_energy = fs->make<TH1F>("num_pu_energy", "num_pu_energy", binsEnergy, 0.f, 300.f);
  num_pu_pt = fs->make<TH1F>("num_pu_pt", "num_pu_pt", binsEnergy, 0.f, 200.f);

  den_sim_eta = fs->make<TH1F>("den_sim_eta", "den_sim_eta", binsEta, -3.0, 3.0);
  den_sim_phi = fs->make<TH1F>("den_sim_phi", "den_sim_phi", binsEta, -3.14, 3.14);
  den_sim_energy = fs->make<TH1F>("den_sim_energy", "den_sim_energy", binsEnergy, 0.f, 300.f);
  den_sim_pt = fs->make<TH1F>("den_sim_pt", "den_sim_pt", binsEnergy, 0.f, 200.f);

  den_reco_eta = fs->make<TH1F>("den_reco_eta", "den_reco_eta", binsEta, -3.0, 3.0);
  den_reco_phi = fs->make<TH1F>("den_reco_phi", "den_reco_phi", binsEta, -3.14, 3.14);
  den_reco_energy = fs->make<TH1F>("den_reco_energy", "den_reco_energy", binsEnergy, 0.f, 300.f);
  den_reco_pt = fs->make<TH1F>("den_reco_pt", "den_reco_pt", binsEnergy, 0.f, 200.f);

  den_fake_reco_eta = fs->make<TH1F>("den_fake_reco_eta", "den_fake_reco_eta", binsEta, -3.0, 3.0);
  den_fake_reco_phi = fs->make<TH1F>("den_fake_reco_phi", "den_fake_reco_phi", binsEta, -3.14, 3.14);
  den_fake_reco_energy = fs->make<TH1F>("den_fake_reco_energy", "den_fake_reco_energy", binsEnergy, 0.f, 300.f);
  den_fake_reco_pt = fs->make<TH1F>("den_fake_reco_pt", "den_fake_reco_pt", binsEnergy, 0.f, 200.f);
  */
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
  desc.add<edm::InputTag>("trackstersMerged", edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("simtrackstersCP", edm::InputTag("ticlSimTracksters", "fromCPs"));
  desc.add<edm::InputTag>("simtrackstersPU", edm::InputTag("ticlSimTracksters", "PU"));
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("recoToSimAssociatorCP",
                          edm::InputTag("tracksterSimTracksterAssociationLinking", "recoToSim"));
  desc.add<edm::InputTag>("simToRecoAssociatorCP",
                          edm::InputTag("tracksterSimTracksterAssociationLinking", "simToReco"));
  desc.add<edm::InputTag>("recoToSimAssociatorPU",
                          edm::InputTag("tracksterSimTracksterAssociationLinking", "recoToSim"));
  descriptions.add("simpleValidation", desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleValidation);
