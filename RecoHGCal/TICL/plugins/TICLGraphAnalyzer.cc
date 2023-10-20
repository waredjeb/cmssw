// Original Authors:  Philipp Zehetner, Wahid Redjeb

#include "TTree.h"
#include "TFile.h"
#include <cassert>

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

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLGraph.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Common.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

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

class TICLGraphAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
public:
  explicit TICLGraphAnalyzer(const edm::ParameterSet&);
  ~TICLGraphAnalyzer() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  typedef math::XYZVector Vector;
  typedef std::vector<double> Vec;

private:
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override;

  void initialize(const HGCalDDDConstants* hgcons,
                  const hgcal::RecHitTools rhtools,
                  const edm::ESHandle<MagneticField> bfieldH,
                  const edm::ESHandle<Propagator> propH);
  void buildLayers();

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override{};
  void endJob() override;

  // Define Tokens
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_token_;
  const edm::EDGetTokenT<TICLGraph> ticlGraph_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> simTracksters_CP_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> simTracksters_PU_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimCP_token_;
  const edm::EDGetTokenT<hgcal::SimToRecoCollectionSimTracksters> tsSimToRecoCP_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimPU_token_;
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloparticles_token_;

  // Output tree
  TTree* tree_;

  void clearVariables();

  TTree* graph_tree;
  std::vector<float> num_eff_energy;
  std::vector<float> num_eff_eta;
  std::vector<float> num_eff_pt;
  std::vector<float> num_eff_phi;
  std::vector<float> den_eff_energy;
  std::vector<float> den_eff_eta;
  std::vector<float> den_eff_pt;
  std::vector<float> den_eff_phi;
  std::vector<float> den_cont_energy;
  std::vector<float> den_cont_eta;
  std::vector<float> den_cont_pt;
  std::vector<float> den_cont_phi;
  std::vector<float> num_cont_energy;
  std::vector<float> num_cont_eta;
  std::vector<float> num_cont_pt;
  std::vector<float> num_cont_phi;
  std::vector<float> fractions;
  int totMerged = 0;
  int totComponents = 0;
  std::vector<int> totNumberOfEdges;

  int ev_event_;
};

void TICLGraphAnalyzer::clearVariables(){
    // event info
};

TICLGraphAnalyzer::TICLGraphAnalyzer(const edm::ParameterSet& ps)
    : tracksters_token_(consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("trackstersclue3d"))),
      ticlGraph_token_(consumes<TICLGraph>(ps.getParameter<edm::InputTag>("ticlGraph"))),
      simTracksters_CP_token_(
          consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("simtrackstersCP"))),
      simTracksters_PU_token_(
          consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("simtrackstersPU"))),
      tsRecoToSimCP_token_(
          consumes<hgcal::RecoToSimCollectionSimTracksters>(ps.getParameter<edm::InputTag>("recoToSimAssociatorCP"))),
      tsSimToRecoCP_token_(
          consumes<hgcal::SimToRecoCollectionSimTracksters>(ps.getParameter<edm::InputTag>("simToRecoAssociatorCP"))),
      tsRecoToSimPU_token_(
          consumes<hgcal::RecoToSimCollectionSimTracksters>(ps.getParameter<edm::InputTag>("recoToSimAssociatorPU"))),
      caloparticles_token_(consumes(ps.getParameter<edm::InputTag>("caloparticles"))){};

TICLGraphAnalyzer::~TICLGraphAnalyzer() { clearVariables(); };

void TICLGraphAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& es) {}

// Define tree and branches
void TICLGraphAnalyzer::beginJob() {
  edm::Service<TFileService> fs;
  graph_tree = fs->make<TTree>("tracksters", "TICL tracksters");
  graph_tree->Branch("event", &ev_event_);

  graph_tree->Branch("num_eff_energy", &num_eff_energy);
  graph_tree->Branch("num_eff_eta", &num_eff_eta);
  graph_tree->Branch("num_eff_pt", &num_eff_pt);
  graph_tree->Branch("num_eff_phi", &num_eff_phi);
  graph_tree->Branch("den_eff_energy", &den_eff_energy);
  graph_tree->Branch("den_eff_eta", &den_eff_eta);
  graph_tree->Branch("den_eff_pt", &den_eff_pt);
  graph_tree->Branch("den_eff_phi", &den_eff_phi);
  graph_tree->Branch("den_cont_energy", &den_cont_energy);
  graph_tree->Branch("den_cont_eta", &den_cont_eta);
  graph_tree->Branch("den_cont_pt", &den_cont_pt);
  graph_tree->Branch("den_cont_phi", &den_cont_phi);
  graph_tree->Branch("num_cont_energy", &num_cont_energy);
  graph_tree->Branch("num_cont_eta", &num_cont_eta);
  graph_tree->Branch("num_cont_pt", &num_cont_pt);
  graph_tree->Branch("num_cont_phi", &num_cont_phi);
  graph_tree->Branch("fractions", &fractions);
  graph_tree->Branch("totNumberOfEdges", &totNumberOfEdges);

  graph_tree->Branch("totMerged", &totMerged);
  graph_tree->Branch("totComponents", &totComponents);
}

void TICLGraphAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  ev_event_ += 1;
  clearVariables();
  //get all the tracksters
  edm::Handle<std::vector<ticl::Trackster>> tracksters_handle;
  event.getByToken(tracksters_token_, tracksters_handle);
  const auto& tracksters = *tracksters_handle;

  edm::Handle<TICLGraph> ticlGraph_handle;
  event.getByToken(ticlGraph_token_, ticlGraph_handle);
  auto& ticlGraph = *ticlGraph_handle;

  //get all the layer clusters
  // simTracksters from CP
  edm::Handle<std::vector<ticl::Trackster>> simTrackstersCP_h;
  event.getByToken(simTracksters_CP_token_, simTrackstersCP_h);
  const auto& simTrackstersCP = *simTrackstersCP_h;

  // simTracksters from PU
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

  edm::Handle<std::vector<CaloParticle>> caloparticles_h;
  event.getByToken(caloparticles_token_, caloparticles_h);
  const auto& caloparticles = *caloparticles_h;

  auto totNumberOfEdgesEv = 0;

  for (auto const& node : ticlGraph.getNodes()) {
    totNumberOfEdgesEv += node.getNeighbours().size();
  }

  totNumberOfEdges.push_back(totNumberOfEdgesEv);

  auto const connectedComponents = ticlGraph.getConnectedComponents();
  int iComp = 0;
  ////std::cout << "Initial graph " << std::endl;
  int count = 0;

  std::vector<int> simTracksterInComponent(connectedComponents.size());
  for (size_t iSim = 0; iSim != simTrackstersCP.size(); iSim++) {
    auto const& simTrackster = simTrackstersCP[iSim];
    count += 1;
    auto const& cp = caloparticles[simTrackster.seedIndex()];
    den_eff_energy.push_back(simTrackster.regressed_energy());
    den_eff_eta.push_back(cp.eta());
    den_eff_phi.push_back(cp.phi());
    den_eff_pt.push_back(cp.pt());

    bool matched = false;
    if (!matched) {
      const edm::Ref<ticl::TracksterCollection> simRef(simTrackstersCP_h, iSim);

      auto const& iSimMap_iter = tsSimToRecoCPMap.find(simRef);

      if (iSimMap_iter != tsSimToRecoCPMap.end()) {
        auto const& iSimMap = iSimMap_iter->val;
        auto iComp = 0;
        for (auto const& comp : connectedComponents) {
          auto sumComponent = 0.;
          auto sumContamination = 0.;
          auto sumEnergy = 0.;
          for (const int iReco : comp) {
            sumEnergy += tracksters[iReco].raw_energy();
            auto const& recoRef = edm::Ref<ticl::TracksterCollection>(tracksters_handle, iReco);
            for (auto const& rts : iSimMap) {
              auto rts_id = (rts.first).get() - (edm::Ref<ticl::TracksterCollection>(tracksters_handle, 0)).get();
              if (iReco == rts_id) {
                sumComponent += rts.second.first;  // shared fraction
                auto const& rtsPU_iter = tsRecoSimPUMap.find(recoRef);
                if (rtsPU_iter != tsRecoSimPUMap.end()) {
                  int iLoop = 0;
                  for (auto const& rtsPU : rtsPU_iter->val) {
                    //                    std::cout << " iLoop " << iLoop << " Trackster energy " << tracksters[iReco].raw_energy()
                    //                              << " rtsPU.second.first " << rtsPU.second.first << std::endl;
                    iLoop += 1;
                    sumContamination += rtsPU.second.first;
                  }
                }
              }
            }
          }
          auto const fraction = sumComponent / simTrackster.raw_energy();
          //         std::cout << "CP energy " << simTrackster.regressed_energy() << " SumEnergy " << sumEnergy << " sumComp "
          //                 << sumComponent << " sumContamination " << sumContamination << std::endl;
          if (sumComponent / simTrackster.raw_energy() >= 0.7 and !matched) {
            num_eff_energy.push_back(simTrackster.regressed_energy());
            num_eff_eta.push_back(cp.eta());
            num_eff_phi.push_back(cp.phi());
            num_eff_pt.push_back(cp.pt());
            fractions.push_back(fraction);
            simTracksterInComponent[iComp] += 1;
            matched = true;
            //std::cout << " Fraction Before " << sumContamination / sumEnergy << " SumCont  " << sumContamination << " sum energy " << sumEnergy << std::endl;
            if (sumContamination / sumEnergy >= 0.4) {
              std::cout << "Sum Contamination " << sumContamination / sumEnergy << std::endl;
              num_cont_energy.push_back(simTrackster.regressed_energy());
              num_cont_eta.push_back(cp.eta());
              num_cont_phi.push_back(cp.phi());
              num_cont_pt.push_back(cp.pt());
            }
          }
          iComp += 1;
        }
      }
    }
  }
  totComponents += simTracksterInComponent.size();
  for (auto const& stsInComp : simTracksterInComponent) {
    if (stsInComp > 1)
      totMerged += 1;
  }
}

void TICLGraphAnalyzer::endJob() { graph_tree->Fill(); }

void TICLGraphAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackstersclue3d", edm::InputTag("ticlTrackstersCLUE3DHigh"));
  desc.add<edm::InputTag>("ticlGraph", edm::InputTag("ticlGraph"));
  desc.add<edm::InputTag>("simtrackstersCP", edm::InputTag("ticlSimTracksters", "fromCPs"));
  desc.add<edm::InputTag>("simtrackstersPU", edm::InputTag("ticlSimTracksters", "PU"));
  desc.add<edm::InputTag>("caloparticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("recoToSimAssociatorCP",
                          edm::InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3D", "recoToSim"));
  desc.add<edm::InputTag>("simToRecoAssociatorCP",
                          edm::InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3D", "simToReco"));
  desc.add<edm::InputTag>("recoToSimAssociatorPU",
                          edm::InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3DPU", "recoToSim"));

  descriptions.add("ticlGraphAnalyzer", desc);
}

DEFINE_FWK_MODULE(TICLGraphAnalyzer);
