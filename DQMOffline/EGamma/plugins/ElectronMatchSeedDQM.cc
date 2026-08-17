/*
Description:
s
DQM module for GSF elecron pixel seeding validation. Measures HLT GSF electron reco efficiency w.r.t. generator-level electrons and monitors the pixel seeds (doublet/triplet counts and seed hit layers). Produces histograms that (should be) viewable in the DQM GUI. 
It fills raw numerator and denominator histograms per event. All divisions with efficiences, fractions, are done at the end via the harvester, electronMatchSeedPostProcessor_cfi.py (a DQMGenericClient), because parallel jobs can sum raw counts but not average ratios. 
Naming conventions for the variables are more detailed but then booked like the standalone analyzer's-
Matching-level histograms:
electronSim(_eta/_phi): all selected gen electrons, gen kinematics (efficiency denominator)
eG_pt/eta/phi: gen kinematics of gen electrons matched by at least 1 reco (efficiency numerator)
eMatched, eR_eta/phi: matched reco electrons, RECO kinematics (kept but wont be used as the efficiency numerator)
electronReco(_eta/_phi): all reconstructed hlt electrons
Seed-level histograms: (port of eDouble/eTrip, initialDoublets/Triplets, BPIX/FPIX):
initialDoublets/initialTriplets/initialQuadPlus: per-event multiplicity of 2-hit/3-hit/>=4-hit seeds
eDouble/eTrip/eQuad: matched electrons (reco pT)
FracFromDoublet/Triplet/QuadPlus_vs_pt = eDouble(eTrip,eQuad)/eMatched
eMatched_BPIX/eMatched_FPIX: layer/disk occupancy of matched-electron seed hits
initialSeedHits: hit multiplicity of every pixel seed (one entry per seed)
eR_seedHits: hit multiplicity of the seed of each matched reco electron
The seed size is not always 2 or 3: seeds built from pixel tracks
(SeedGeneratorFromProtoTracksEDProducer with includeFourthHit=True, as enabled by the
hltEgammaPixelTrackSeeding process modifier) routinely carry 4 or more hits. The
initialSeedHits/eR_seedHits distributions cover the full range, so the doublet/triplet/
quad-plus categories always add up to the total.
*/
#include <memory>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
//Sim-truth (gen-level) pixel occupancy
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include <set>
#include <map>
//In CMSSW
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
//Seed-level added includes
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
namespace {
  //seed sizes are plotted with integer bin centres 0..10; anything larger is piled
  //into the last bin rather than lost to the overflow
  constexpr unsigned int kMaxSeedNHitsBin = 10;
  double clampSeedNHits(unsigned int nHits) { return std::min(nHits, kMaxSeedNHitsBin); }
}  // namespace
class ElectronMatchSeedDQM : public DQMEDAnalyzer {
public:
  explicit ElectronMatchSeedDQM(const edm::ParameterSet&);
  ~ElectronMatchSeedDQM() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
//private:
//book histograms is called oncce per run to create the MonitorElements
//analyze is called once per event to fill
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  void analyze(edm::Event const&, edm::EventSetup const&) override;
private:
  //helper methods for matching logic (from ElectronMatchSeedNew)
  reco::GenParticle get_lastcopy_prefsr(reco::GenParticle part);
  reco::GenParticleCollection get_genparts(const reco::GenParticleCollection& genparts);
  //genparticles changed to selected below
  reco::GenParticle const* match_to_gen(double eta_reco, double phi_reco, const reco::GenParticleCollection& selected);
  //tokens
  const edm::EDGetTokenT<reco::GenParticleCollection> gensToken_;
  const edm::EDGetTokenT<reco::ElectronCollection> electronCollectionToken_;
  const edm::EDGetTokenT<reco::ElectronSeedCollection> pixelSeedsToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const edm::EDGetTokenT<edm::SimTrackContainer> simTracksToken_;
  std::vector<edm::EDGetTokenT<edm::PSimHitContainer>> pixelSimHitTokens_;
  //parameters
  double maxDeltaR_;
  std::string folderName_;
  //matching-level MonitorElements
  MonitorElement* h_electronSim_; //denominator: total gen electrons pT
  MonitorElement* h_electronSimEta_;
  MonitorElement* h_electronSimPhi_;
  MonitorElement* h_electronReco_; //total reconstructed electrons pT
  MonitorElement* h_electronRecoEta_;
  MonitorElement* h_electronRecoPhi_;
  MonitorElement* h_eMatched_; 
  MonitorElement* h_eMatchedPt_;  
  MonitorElement* h_eMatchedEta_;
  MonitorElement* h_eMatchedPhi_;
  MonitorElement* h_eMatchedGenPt_; //numerator: gen kinematics of the match (eG_* in tree)
  MonitorElement* h_eMatchedGenEta_;
  MonitorElement* h_eMatchedGenPhi_;
  //seed-level MonitorElements
  MonitorElement* h_nSeedsTotal_; //per-event total seed count
  MonitorElement* h_nSeedsDoublets_; //per-event count of 2-hit seeds (initialDoublets)
  MonitorElement* h_nSeedsTriplets_; //per-event count of 3-hit seeds (initialTriplets)
  MonitorElement* h_nSeedsQuadPlus_; //per-event count of >=4-hit seeds (initialQuadPlus)
  MonitorElement* h_seedNHits_; //hit multiplicity of every seed, one entry per seed (initialSeedHits)
  MonitorElement* h_eMatchedFromDoublet_; //matched reco pt, seed had 2 hits (eDouble)
  MonitorElement* h_eMatchedFromTriplet_; //matched reco pt, seed had 3 hits (eTrip)
  MonitorElement* h_eMatchedFromQuadPlus_; //matched reco pt, seed had >=4 hits (eQuad)
  MonitorElement* h_eMatchedSeedNHits_; //hit multiplicity of the matched electron's seed (eR_seedHits)
  MonitorElement* h_seedHitsBPIXLayer_; //barrel pixel layer of matched-electron seed hits (eMatched_BPIX)
  MonitorElement* h_seedHitsFPIXDisk_; //forward pixel disk of matched-electron seed hits (eMatched_FPIX)
  MonitorElement* h_eGenBPIXLayer_; //BPIX layers crossed by gen electrons (SimHits)
  MonitorElement* h_eGenFPIXDisk_; //FPIX disks crossed by gen electrons (SimHits)
  //total seeds vs. All Gen Electrons (Denominator)
  MonitorElement* h_nSeedsVsGenPt_;
  MonitorElement* h_nSeedsVsGenEta_;
  MonitorElement* h_nSeedsVsGenPhi_;
  //total seeds vs. Reco Electrons
  MonitorElement* h_nSeedsVsRecoPt_;
  MonitorElement* h_nSeedsVsRecoEta_;
  MonitorElement* h_nSeedsVsRecoPhi_;
  //total seeds vs. Matched Gen Electrons (Numerator)
  MonitorElement* h_nSeedsVsMatchedGenPt_;
  MonitorElement* h_nSeedsVsMatchedGenEta_;
  MonitorElement* h_nSeedsVsMatchedGenPhi_;
};
ElectronMatchSeedDQM::ElectronMatchSeedDQM(const edm::ParameterSet& iConfig)
    : gensToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genP"))),
      electronCollectionToken_(consumes<reco::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
      pixelSeedsToken_(consumes<reco::ElectronSeedCollection>(iConfig.getParameter<edm::InputTag>("pixelSeedsProducer"))),
      topoToken_(esConsumes()),
      simTracksToken_(consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("simTracks"))),
      maxDeltaR_(iConfig.getParameter<double>("DeltaR")),
      folderName_(iConfig.getParameter<std::string>("dqmFolder")) {
    for (const auto& tag : iConfig.getParameter<std::vector<edm::InputTag>>("pixelSimHits"))
        pixelSimHitTokens_.push_back(consumes<edm::PSimHitContainer>(tag));
      }
void ElectronMatchSeedDQM::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const& iRun, edm::EventSetup const& iSetup) {
  //set top-level DQM subsystem directory path
  ibooker.setCurrentFolder("HLT/EGamma/" + folderName_);
  //pt: 30 bins 0-150 (matches standalone)
  //numerator/denominator must share binning
  //note: numerator and denominator must share binning for the division in the harvester to work
  h_electronSim_ = ibooker.book1D("electronSim", "gen electron p_{T};p_{T} [GeV];entries", 30, 0, 150);
  h_electronReco_ = ibooker.book1D("electronReco", "reco electron p_{T};p_{T} [GeV];entries", 30, 0, 150);
  h_eMatched_ = ibooker.book1D("eMatched", "matched reco electron p_{T};p_{T} [GeV];entries", 30, 0, 150);
  //eta/phi: same variables for eff vs. eta/phi in the harvester
  h_electronSimEta_ = ibooker.book1D("electronSim_eta", "gen electron #eta;#eta;entries", 30, -3.0, 3.0);
  h_electronSimPhi_ = ibooker.book1D("electronSim_phi", "gen electron #phi;#phi;entries", 30, -3.1416, 3.1416);
  h_electronRecoEta_ = ibooker.book1D("electronReco_eta", "reco electron #eta;#eta;entries", 30, -3.0, 3.0);
  h_electronRecoPhi_ = ibooker.book1D("electronReco_phi", "reco electron #phi;#phi;entries", 30, -3.1416, 3.1416);
  h_eMatchedPt_ = ibooker.book1D("eR_pt", "matched reco electron p_{T};p_{T} [GeV];entries", 30, 0, 150);
  h_eMatchedEta_ = ibooker.book1D("eR_eta", "matched reco electron #eta;#eta;entries", 30, -3.0, 3.0);
  h_eMatchedPhi_ = ibooker.book1D("eR_phi", "matched reco electron #phi;#phi;entries", 30, -3.1416, 3.1416);
  //gen kinematics of the matched partner (tree branches eG_pt/eta/phi)
  h_eMatchedGenPt_ = ibooker.book1D("eG_pt", "gen p_{T};p_{T} [GeV];entries", 30, 0, 150);
  h_eMatchedGenEta_ = ibooker.book1D("eG_eta", "gen #eta of matched electron;#eta;entries", 30, -3.0, 3.0);
  h_eMatchedGenPhi_ = ibooker.book1D("eG_phi", "gen #phi of matched electron;#phi;entries", 30, -3.1416, 3.1416);
  //seed-level (binning confirm ranges with seeding group)
  h_nSeedsTotal_ = ibooker.book1D("initialSeeds", "pixel seeds per event;n seeds;events", 100, 0, 100); 
  h_nSeedsDoublets_ = ibooker.book1D("initialDoublets", "pixel doublet seeds per event;n doublet seeds;events", 50, 0, 50);
  h_nSeedsTriplets_ = ibooker.book1D("initialTriplets", "pixel triplet seeds per event;n triplet seeds;events", 50, 0, 50);
  h_nSeedsQuadPlus_ = ibooker.book1D("initialQuadPlus", "pixel seeds with #geq4 hits per event;n #geq4-hit seeds;events", 50, 0, 50);
  h_eMatchedFromDoublet_ = ibooker.book1D("eDouble", "matched electrons from doublet seeds;p_{T} [GeV];entries", 30, 0, 150);
  h_eMatchedFromTriplet_ = ibooker.book1D("eTrip", "matched electrons from triplet seeds;p_{T} [GeV];entries", 30, 0, 150);
  h_eMatchedFromQuadPlus_ = ibooker.book1D("eQuad", "matched electrons from #geq4-hit seeds;p_{T} [GeV];entries", 30, 0, 150);
  //seed size: one entry per seed / per matched electron. Integer bin centres 0..10,
  //the last bin collects any seed with more hits than that.
  h_seedNHits_ = ibooker.book1D("initialSeedHits", "pixel seed hit multiplicity;n hits on seed;seeds", 11, -0.5, 10.5);
  h_eMatchedSeedNHits_ = ibooker.book1D("eR_seedHits", "matched-electron seed hit multiplicity;n hits on seed;electrons", 11, -0.5, 10.5);
  //seed hits of matched reco electrons: one entry per HIT
  h_seedHitsBPIXLayer_ = ibooker.book1D("eR_BPIX", "matched-electron seed hits: BPIX layer;layer;hits", 4, 0.5, 4.5);
  h_seedHitsFPIXDisk_ = ibooker.book1D("eR_FPIX", "matched-electron seed hits: FPIX disk;disk;hits", 12, 0.5, 12.5);
  //total seeds vs. gen kinematics
  h_nSeedsVsGenPt_  = ibooker.bookProfile("nSeeds_vs_genPt",  "Total seeds vs gen p_{T};gen p_{T} [GeV];n seeds", 30, 0, 150, 50, 0, 50);
  h_nSeedsVsGenEta_ = ibooker.bookProfile("nSeeds_vs_genEta", "Total seeds vs gen #eta;gen #eta;n seeds",       30, -3.0, 3.0, 50, 0, 50);
  h_nSeedsVsGenPhi_ = ibooker.bookProfile("nSeeds_vs_genPhi", "Total seeds vs gen #phi;gen #phi;n seeds",       30, -3.1416, 3.1416, 50, 0, 50);
  //total seeds vs. reco kinematics
  h_nSeedsVsRecoPt_  = ibooker.bookProfile("nSeeds_vs_recoPt",  "Total seeds vs reco p_{T};reco p_{T} [GeV];n seeds", 30, 0, 150, 50, 0, 50);
  h_nSeedsVsRecoEta_ = ibooker.bookProfile("nSeeds_vs_recoEta", "Total seeds vs reco #eta;reco #eta;n seeds",       30, -3.0, 3.0, 50, 0, 50);
  h_nSeedsVsRecoPhi_ = ibooker.bookProfile("nSeeds_vs_recoPhi", "Total seeds vs reco #phi;reco #phi;n seeds",       30, -3.1416, 3.1416, 50, 0, 50);
  //total seeds vs. matched gen kinematics
  h_nSeedsVsMatchedGenPt_  = ibooker.bookProfile("nSeeds_vs_matchedGenPt",  "Total seeds vs matched gen p_{T};gen p_{T} [GeV];n seeds", 30, 0, 150, 50, 0, 50);
  h_nSeedsVsMatchedGenEta_ = ibooker.bookProfile("nSeeds_vs_matchedGenEta", "Total seeds vs matched gen #eta;gen #eta;n seeds",       30, -3.0, 3.0, 50, 0, 50);
  h_nSeedsVsMatchedGenPhi_ = ibooker.bookProfile("nSeeds_vs_matchedGenPhi", "Total seeds vs matched gen #phi;gen #phi;n seeds",       30, -3.1416, 3.1416, 50, 0, 50);
  //sim-truth occupancy: which pixel layers the selected gen electrons actually crossed,
  //Deduplicated per electron, so one entry per electron per layer.
  //Same binning as eR_BPIX/eR_FPIX 
  h_eGenBPIXLayer_ = ibooker.book1D("eG_BPIX", "gen electron sim hits: BPIX layer;layer;electrons", 4, 0.5, 4.5);
  h_eGenFPIXDisk_ = ibooker.book1D("eG_FPIX", "gen electron sim hits: FPIX disk;disk;electrons", 12, 0.5, 12.5);
}
void ElectronMatchSeedDQM::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  auto genParticleH = iEvent.getHandle(gensToken_);
  if (!genParticleH.isValid())
    return;
    const reco::GenParticleCollection& genParticles = *genParticleH;
    auto simEles = get_genparts(genParticles);
    auto pixelSeedsH = iEvent.getHandle(pixelSeedsToken_);
    int nTotal = 0;
    int nDoublets = 0;
    int nTriplets = 0;
    int nQuadPlus = 0;
    if (pixelSeedsH.isValid()) {
      nTotal = pixelSeedsH->size();
      for (const auto& seed : *pixelSeedsH) {
        h_seedNHits_->Fill(clampSeedNHits(seed.nHits()));
        if (seed.nHits() == 2)
          ++nDoublets;
        else if (seed.nHits() == 3)
          ++nTriplets;
        else if (seed.nHits() >= 4)
          ++nQuadPlus;
      }
    }
    h_nSeedsTotal_->Fill(nTotal);
    h_nSeedsDoublets_->Fill(nDoublets);
    h_nSeedsTriplets_->Fill(nTriplets);
    h_nSeedsQuadPlus_->Fill(nQuadPlus);
  //fill gen electron histograms
  for (const auto& simEle : simEles) {
    h_electronSim_->Fill(simEle.pt());
    h_electronSimEta_->Fill(simEle.eta());
    h_electronSimPhi_->Fill(simEle.phi());
    //fill total seeds vs. all gen kinematics
    h_nSeedsVsGenPt_->Fill(simEle.pt(), nTotal);
    h_nSeedsVsGenEta_->Fill(simEle.eta(), nTotal);
    h_nSeedsVsGenPhi_->Fill(simEle.phi(), nTotal);
  }
  //TrackerTopology is needed by both the gen-level and reco-level blocks below
  const TrackerTopology& ttopo = iSetup.getData(topoToken_);
  //gen-level pixel occupancy from GEANT4 truth. SimTracks are dR-matched to the selected
  //gen electrons, then every pixel SimHit carrying a matched trackId contributes its
  //layer/disk. Deduplicated per electron: crossing layer 1 twice counts once.
  auto simTracksH = iEvent.getHandle(simTracksToken_);
  if (simTracksH.isValid()) {
    //trackId -> index of the gen electron it belongs to
    std::map<unsigned int, size_t> trackIdToGen;
    for (const auto& st : *simTracksH) {
      if (std::abs(st.type()) != 11)
        continue;
      for (size_t g = 0; g < simEles.size(); ++g) {
        double dr2 = reco::deltaR2(st.momentum().eta(), st.momentum().phi(), simEles[g].eta(), simEles[g].phi());
        if (dr2 < maxDeltaR_ * maxDeltaR_) {
          trackIdToGen[st.trackId()] = g;
          break;
        }
      }
    }
    std::vector<std::set<int>> genBpix(simEles.size());
    std::vector<std::set<int>> genFpix(simEles.size());
    for (const auto& tok : pixelSimHitTokens_) {
      auto hitsH = iEvent.getHandle(tok);
      if (!hitsH.isValid())
        continue;
      for (const auto& hit : *hitsH) {
        auto it = trackIdToGen.find(hit.trackId());
        if (it == trackIdToGen.end())
          continue;
        DetId det(hit.detUnitId());
        if (det.det() != DetId::Tracker)
          continue;
        int subdet = det.subdetId();
        if (subdet == PixelSubdetector::PixelBarrel)
          genBpix[it->second].insert(ttopo.pxbLayer(det));
        else if (subdet == PixelSubdetector::PixelEndcap)
          genFpix[it->second].insert(ttopo.pxfDisk(det));
      }
    }
    for (size_t g = 0; g < simEles.size(); ++g) {
      for (int l : genBpix[g])
        h_eGenBPIXLayer_->Fill(l);
      for (int d : genFpix[g])
        h_eGenFPIXDisk_->Fill(d);
    }
  }
  auto gsfElectronH = iEvent.getHandle(electronCollectionToken_);
    if (gsfElectronH.isValid()) {
      std::vector<bool> genMatched(simEles.size(), false);
      for (const auto& gsfElectron : *gsfElectronH) {
        //all reco electrons: electronReco + eR_pt (fills eR_pt for every reco)
        h_electronReco_->Fill(gsfElectron.pt());
        h_electronRecoEta_->Fill(gsfElectron.eta());
        h_electronRecoPhi_->Fill(gsfElectron.phi());
        h_eMatchedPt_->Fill(gsfElectron.pt());
        //fill total seeds vs. reco kinematics
        h_nSeedsVsRecoPt_->Fill(gsfElectron.pt(), nTotal);
        h_nSeedsVsRecoEta_->Fill(gsfElectron.eta(), nTotal);
        h_nSeedsVsRecoPhi_->Fill(gsfElectron.phi(), nTotal);
        auto matchedEle = match_to_gen(gsfElectron.eta(), gsfElectron.phi(), simEles);
        if (matchedEle == nullptr)
          continue;   //unmatched electrons stop here
        genMatched[matchedEle - simEles.data()] = true;
        //MATCHED reco electrons only
        h_eMatched_->Fill(gsfElectron.pt());
        h_eMatchedEta_->Fill(gsfElectron.eta());
        h_eMatchedPhi_->Fill(gsfElectron.phi());
        //seed classification (eDouble/eTrip, BPIX/FPIX) 
        auto gsfTrack = gsfElectron.gsfTrack();
        if (gsfTrack.isNull())
          continue;
        auto seedRef = gsfTrack->seedRef();
        if (seedRef.isNull() || !seedRef.isAvailable())
          continue;
        unsigned int nHits = seedRef->nHits();
        h_eMatchedSeedNHits_->Fill(clampSeedNHits(nHits));
        if (nHits == 2)
          h_eMatchedFromDoublet_->Fill(gsfElectron.pt());
        else if (nHits == 3)
          h_eMatchedFromTriplet_->Fill(gsfElectron.pt());
        else if (nHits >= 4)
          h_eMatchedFromQuadPlus_->Fill(gsfElectron.pt());
        for (const auto& rhit : seedRef->recHits()) {
          if (!rhit.isValid())
            continue;
          DetId det = rhit.geographicalId();
          if (det.det() != DetId::Tracker)
            continue;
          int subdet = det.subdetId();
          if (subdet == PixelSubdetector::PixelBarrel)
            h_seedHitsBPIXLayer_->Fill(ttopo.pxbLayer(det));
          else if (subdet == PixelSubdetector::PixelEndcap)
            h_seedHitsFPIXDisk_->Fill(ttopo.pxfDisk(det));
        }
      }
      //efficiency numerator: once per gen electron matched by at least one reco
      for (size_t i = 0; i < simEles.size(); ++i) {
        if (genMatched[i]) {
          h_eMatchedGenPt_->Fill(simEles[i].pt());
          h_eMatchedGenEta_->Fill(simEles[i].eta());
          h_eMatchedGenPhi_->Fill(simEles[i].phi());
          //fill total seeds vs. matched gen kinematics
          h_nSeedsVsMatchedGenPt_->Fill(simEles[i].pt(), nTotal);
          h_nSeedsVsMatchedGenEta_->Fill(simEles[i].eta(), nTotal);
          h_nSeedsVsMatchedGenPhi_->Fill(simEles[i].phi(), nTotal);
        }
      }
    }
  }
//ALL changes are from this point and above for seed-level
reco::GenParticle ElectronMatchSeedDQM::get_lastcopy_prefsr(reco::GenParticle part) {
  auto daughters = part.daughterRefVector();
  if (daughters.size() == 1 && daughters.at(0)->pdgId() == part.pdgId()) {
    return get_lastcopy_prefsr(*(daughters.at(0)));
  }
  return part;
}
reco::GenParticleCollection ElectronMatchSeedDQM::get_genparts(const reco::GenParticleCollection& genparts) {
  std::vector<reco::GenParticle> selected;
  for (const auto& part : genparts) {
    auto pdg_id = part.pdgId();
    if (abs(pdg_id) == 11) {
      if (part.numberOfMothers() == 0) {
        selected.push_back(part);
        continue;
      }
      if (part.mother() != nullptr && part.isHardProcess()) {
        auto momId = fabs(part.mother()->pdgId());
        if (momId == 24 || momId == 23) {
          selected.push_back(get_lastcopy_prefsr(part));
        }
      }
    }
  }
  return selected;
}
reco::GenParticle const* ElectronMatchSeedDQM::match_to_gen(double eta_reco,
                                                            double phi_reco,
                                                            const reco::GenParticleCollection& selected) {
  reco::GenParticle const* best_match = nullptr;
  double best_dr2 = maxDeltaR_ * maxDeltaR_;
  for (const auto& p : selected) {
    double dr2 = reco::deltaR2(eta_reco, phi_reco, p.eta(), p.phi());
    if (dr2 < best_dr2) {
      best_dr2 = dr2;
      best_match = &p;
    }
  }
  return best_match;
}
//below is to check if they match
/*
if (best_match != nullptr) {
    std::cout << "(dqm) match: reco(eta,phi)=(" << eta_reco << "," << phi_reco
              << ") dR=" << std::sqrt(best_dr2) << " (cut=" << maxDeltaR_ << ")" << std::endl;
  } else {
    std::cout << "(dqm) no match: reco(eta,phi)=(" << eta_reco << "," << phi_reco
              << ") cut=" << maxDeltaR_ << std::endl;
  }
  return best_match;
}
*/
void ElectronMatchSeedDQM::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("genP", edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("electrons", edm::InputTag("hltEgammaGsfElectronsL1Seeded"));
  desc.add<edm::InputTag>("pixelSeedsProducer", edm::InputTag("hltEgammaElectronPixelSeedsL1Seeded"));
  desc.add<edm::InputTag>("simTracks", edm::InputTag("g4SimHits"));
  desc.add<std::vector<edm::InputTag>>("pixelSimHits",
                                       {edm::InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"),
                                        edm::InputTag("g4SimHits", "TrackerHitsPixelBarrelHighTof"),
                                        edm::InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof"),
                                        edm::InputTag("g4SimHits", "TrackerHitsPixelEndcapHighTof")});
  desc.add<double>("DeltaR", 0.1);
  desc.add<std::string>("dqmFolder", "PixelSeedMatching");
  descriptions.add("electronMatchSeedDQM", desc);
}
DEFINE_FWK_MODULE(ElectronMatchSeedDQM);