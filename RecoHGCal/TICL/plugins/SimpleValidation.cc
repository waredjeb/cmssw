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
  int global_reco_fake_ = 0;
  int global_reco_fake_high_energy_ = 0;
  int global_reco_fake_low_energy_ = 0;
  int global_fake_ = 0;
  int global_fake_high_energy_num_ = 0;
  int global_fake_low_energy_num_ = 0;
  int global_fake_high_energy_den_ = 0;
  int global_fake_low_energy_den_ = 0;
  int global_ass_pur_ = 0;
  int global_ass_eff_ = 0;
  int global_merge_ = 0;
  int global_sim_ = 0;
  int pu_reco_ = 0;
  float reco_sim_ratio = 0.f;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_token_;
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloParticles_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> simTracksters_CP_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimCP_token_;
  const edm::EDGetTokenT<hgcal::SimToRecoCollectionSimTracksters> tsSimToRecoCP_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimPU_token_;
  //  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimPU_token_;

  TTree* output_tree_;
};

SimpleValidation::SimpleValidation(const edm::ParameterSet& iConfig)
    : tracksters_token_(
          consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("trackstersclue3d"))),
      caloParticles_token_(consumes<std::vector<CaloParticle>>(iConfig.getParameter<edm::InputTag>("caloParticles"))),
      simTracksters_CP_token_(
          consumes<std::vector<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("simtrackstersCP"))),
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

  // trackster reco to sim CP
  edm::Handle<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimCP_h;
  event.getByToken(tsRecoToSimCP_token_, tsRecoToSimCP_h);
  auto const& tsRecoSimCPMap = *tsRecoToSimCP_h;

  // sim simTrackster CP to reco trackster
  edm::Handle<hgcal::SimToRecoCollectionSimTracksters> tsSimToRecoCP_h;
  event.getByToken(tsSimToRecoCP_token_, tsSimToRecoCP_h);
  auto const& tsSimToRecoCPMap = *tsSimToRecoCP_h;

  edm::Handle<std::vector<CaloParticle>> caloParticles_h;
  event.getByToken(caloParticles_token_, caloParticles_h);
  auto const& caloParticles = *caloParticles_h;

  reco_sim_ratio = tracksters.size() / caloParticles.size();

  edm::Handle<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimPU_h;
  event.getByToken(tsRecoToSimPU_token_, tsRecoToSimPU_h);
  auto const& tsRecoSimPUMap = *tsRecoToSimPU_h;

  std::vector<int> stsInTrackster(tracksters.size(), -1);
  std::vector<int> stsInTracksterSignal(tracksters.size(), 0);
  std::vector<size_t> cPIndices;
  //  std::cout << "caloParticles size " << caloParticles.size() << std::endl;
  removeCPFromPU(caloParticles, cPIndices, true);
  //  std::cout << "cpIndices " << cPIndices.size() << std::endl;

  for (size_t iReco = 0; iReco != tracksters.size(); iReco++) {
    const edm::Ref<ticl::TracksterCollection> tsRef(tracksters_handle, iReco);
    auto const sts_iter = tsRecoSimCPMap.find(tsRef);
    auto const sts_iterPU = tsRecoSimPUMap.find(tsRef);
    if (sts_iterPU != tsRecoSimPUMap.end()) {
      const auto& stsPUAssociated = sts_iterPU->val;
      if (stsPUAssociated[0].second.first / tracksters[iReco].raw_energy() <= 0.95) {
        global_reco_fake_ += 1;
        if (stsInTrackster[iReco] == -1) {
          stsInTrackster[iReco] = 0;
        }
        if (sts_iter != tsRecoSimCPMap.end()) {
          const auto& stsAssociated = sts_iter->val;
          for (auto const& sts : stsAssociated) {
            //auto sts_idx = (sts.first).get() - (edm::Ref<ticl::TracksterCollection>(simTrackstersCP_h, 0)).get();
            if (sts.second.second <= 0.6) {
              stsInTrackster[iReco] += 1;
            }
          }
        }
      }
    }
  }

  for (size_t iReco = 0; iReco != stsInTrackster.size(); iReco++) {
    if (stsInTrackster[iReco] == -1) {
      continue;
    }
    if (stsInTrackster[iReco] > 1) {
      global_merge_ += 1;
    }
    if (tracksters[iReco].raw_energy() >= 25.f) {
      global_fake_high_energy_den_ += 1;
      if (stsInTrackster[iReco] == 0) {
        global_fake_high_energy_num_ += 1;
      }
    } else {
      global_fake_low_energy_den_ += 1;
      if (stsInTrackster[iReco] == 0) {
        global_fake_low_energy_num_ += 1;
      }
    }
  }

  for (size_t iSim = 0; iSim != simTrackstersCP.size(); iSim++) {
    bool matchedPur = false;
    bool matchedEff = false;
    int totMatched = 0;
    bool merged = 0;
    auto cpIndex = simTrackstersCP[iSim].seedIndex();
    if (std::find(cPIndices.begin(), cPIndices.end(), cpIndex) == cPIndices.end()) {
      continue;
    }
    const edm::Ref<ticl::TracksterCollection> stsCPRef(simTrackstersCP_h, iSim);
    auto const ts_iter = tsSimToRecoCPMap.find(stsCPRef);
    if (ts_iter != tsSimToRecoCPMap.end()) {
      const auto& tsAssociated = ts_iter->val;
      for (auto const& ts : tsAssociated) {
        auto ts_idx = (ts.first).get() - (edm::Ref<ticl::TracksterCollection>(tracksters_handle, 0)).get();
        auto const& recoRef = edm::Ref<ticl::TracksterCollection>(tracksters_handle, ts_idx);
        if (ts.second.second <= 0.2 and !matchedPur) {
          global_ass_pur_ += 1;
          matchedPur = true;
        }
        if (ts.second.first / simTrackstersCP[iSim].raw_energy() >= 0.5 and !matchedEff) {
          // check PU contamination
          auto const& pu_iter = tsRecoSimPUMap.find(recoRef);
          if (pu_iter != tsRecoSimPUMap.end()) {
            auto const& puPair = pu_iter->val[0];
            if (puPair.second.first / recoRef->raw_energy() <= 0.2) {
              global_ass_eff_ += 1;
              matchedEff = true;
            }
          }
        }
      }
    }
  }
  global_reco_ += tracksters.size();
  global_sim_ += simTrackstersCP.size();
}

// ------------ method called once each job just before starting event loop  ------------
void SimpleValidation::beginJob() {
  // please remove this method if not needed
  edm::Service<TFileService> fs;
  output_tree_ = fs->make<TTree>("output", "putput params");

  output_tree_->Branch("reco", &global_reco_);
  output_tree_->Branch("recoForFake", &global_reco_fake_);
  output_tree_->Branch("recoForFakeHE", &global_reco_fake_high_energy_);
  output_tree_->Branch("recoForFakeLE", &global_reco_fake_low_energy_);
  output_tree_->Branch("recoAssPur", &global_ass_pur_);
  output_tree_->Branch("reco_ass_eff_", &global_ass_eff_);
  output_tree_->Branch("merge", &global_merge_);
  output_tree_->Branch("pu", &pu_reco_);
  output_tree_->Branch("fake", &global_fake_);
  output_tree_->Branch("fakeHEDen", &global_fake_high_energy_den_);
  output_tree_->Branch("fakeLEDen", &global_fake_low_energy_den_);
  output_tree_->Branch("fakeHENum", &global_fake_high_energy_num_);
  output_tree_->Branch("fakeLENum", &global_fake_low_energy_num_);
  output_tree_->Branch("sim", &global_sim_);
  output_tree_->Branch("reco_sim_ratio", &reco_sim_ratio);
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
  desc.add<edm::InputTag>("trackstersclue3d", edm::InputTag("mergedTracksters"));
  desc.add<edm::InputTag>("simtrackstersCP", edm::InputTag("ticlSimTracksters", "fromCPs"));
  desc.add<edm::InputTag>("caloParticles", edm::InputTag("mix", "MergedCaloTruth"));
  //  desc.add<edm::InputTag>("simtrackstersPU", edm::InputTag("ticlSimTracksters", "PU"));
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("recoToSimAssociatorCP",
                          edm::InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3D", "recoToSim"));
  desc.add<edm::InputTag>("simToRecoAssociatorCP",
                          edm::InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3D", "simToReco"));
  desc.add<edm::InputTag>("recoToSimAssociatorPU",
                          edm::InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3DPU", "recoToSim"));
  descriptions.add("simpleValidation", desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleValidation);
