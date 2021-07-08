#include <numeric>
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "SimDataFormats/Associations/interface/LayerClusterToSimClusterAssociator.h"
#include "SimDataFormats/Associations/interface/LayerClusterToCaloParticleAssociator.h"
#include "SimDataFormats/Associations/interface/TracksterToSimTracksterAssociator.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Vector.h"
#include "TTree.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Dataformats
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/RectangularEtaPhiRegion.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace ticl;

class TrackstersLinkingDebugger : public edm::EDAnalyzer {
public:
  TrackstersLinkingDebugger(const edm::ParameterSet&);
  ~TrackstersLinkingDebugger(){};
  void beginRun(const edm::Run&, const edm::EventSetup&) override;
  void beginJob();
  void endJob();
  void analyze(const edm::Event& event, const edm::EventSetup& es);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void endjob();
  bool isCP(Trackster& t, CaloParticle& cp, edm::ProductID cp_id, edm::Handle<std::vector<CaloParticle>> cp_handle);
  void tracksterInfo(const Trackster& t, const edm::Event& event) const;
  std::pair<std::vector<Trackster>, std::vector<int>> findCompatibleTrackster(std::vector<Trackster> tracksters,
                                                                              const double deltaR,
                                                                              const double min_cos,
                                                                              int tShift,
                                                                              int Nsize,
                                                                              std::vector<int>& mask);
  bool isNuclearInteraction(
      std::vector<Trackster> tracksters, double dR, int t_id1, int trackstersAvailable, std::vector<int>& mask);
  const int associateTrack(std::vector<Trackster> compatibleTrackster, std::vector<reco::Track> tracks);
  math::XYZTLorentzVector assignP4(Trackster t, std::vector<reco::Track> tracks);

  struct Candidate {
    std::vector<Trackster> linkedTrackster;
    std::vector<int> linkedTracksterId;
    reco::Track matchedTracks;
    int matchedTracksId = -99;
  };

  void clear();

private:
  const edm::EDGetTokenT<std::vector<Trackster>> tracksters_token_;
  const edm::EDGetTokenT<std::vector<Trackster>> simtracksters_token_;
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> layerClusters_token;
  const edm::EDGetTokenT<reco::GenParticleCollection> gensToken_;  //GenParticles
  const edm::EDGetTokenT<std::unordered_map<DetId, const HGCRecHit*>> hitmap_token_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometry_token_;
  const edm::EDGetTokenT<std::vector<reco::Track>> tracks_token_;
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloparticles_token_;
  const edm::EDGetTokenT<std::vector<SimCluster>> clusters_token_;

  const double min_angle_;
  const bool do2D;
  TTree* tree_;
  hgcal::RecHitTools rhtools_;

  // std::vector<double>* best_dot = new std::vector<double>;
  // std::vector<double>* genEnergy = new std::vector<double>;
  // std::vector<double>* genEnergySim = new std::vector<double>;
  // std::vector<double>* tracksterSumEnergy = new std::vector<double>;
  // std::vector<double>* tracksterSumEnergyRegressed = new std::vector<double>;
  // std::vector<double>* tracksterSimSumEnergy = new std::vector<double>;
  // std::vector<double>* tracksterSimSumEnergyRegressed = new std::vector<double>;
  // std::vector<double>* simTrackstersEnergy = new std::vector<double>;
  // std::vector<double>* simTrackstersRawEnergy = new std::vector<double>;

  std::vector<double>* cp_pt = new std::vector<double>;
  std::vector<double>* cp_E = new std::vector<double>;
  std::vector<double>* simT_pt = new std::vector<double>;
  std::vector<double>* simT_E = new std::vector<double>;
  std::vector<double>* candidate_pt = new std::vector<double>;
  std::vector<double>* candidate_E = new std::vector<double>;
  std::vector<double>* gen_pt = new std::vector<double>;
  std::vector<double>* gen_E = new std::vector<double>;
  std::vector<bool>* fullMatch = new std::vector<bool>;


  double total_cp = 0.;
  double tot_linkable_trackster = 0.;
  double tot_trackster_linked = 0.;
  double tot_missed_link = 0.;
  double real_double = 0.;
  double tot_alone = 0.;
  double evt_count = 0.;
  double zeroT = 0.;
  double zero = 0;
  double tot_real_alone = 0.;
  double tot_cand_found = 0.;
};

TrackstersLinkingDebugger::TrackstersLinkingDebugger(const edm::ParameterSet& ps)
    : tracksters_token_(consumes<std::vector<Trackster>>(ps.getParameter<edm::InputTag>("tracksters"))),
      simtracksters_token_(consumes<std::vector<Trackster>>(ps.getParameter<edm::InputTag>("simTracksters"))),
      layerClusters_token(consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layerClusters"))),
      gensToken_(consumes<reco::GenParticleCollection>(ps.getParameter<edm::InputTag>("genP"))),
      // hitmap_token_(consumes<std::unordered_map<DetId,const HGCRecHit*>>(ps.getParameter<edm::InputTag>("hitMap"))),
      caloGeometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>()),
      tracks_token_(consumes<std::vector<reco::Track>>(ps.getParameter<edm::InputTag>("tracks"))),
      caloparticles_token_(consumes<std::vector<CaloParticle>>(ps.getParameter<edm::InputTag>("caloparticles"))),
      clusters_token_(consumes<std::vector<SimCluster>>(ps.getParameter<edm::InputTag>("simCluster"))),
      min_angle_(ps.getParameter<double>("min_angle")),
      do2D(ps.getParameter<bool>("do2Ddot")) {}

void TrackstersLinkingDebugger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("tracksters", edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("simTracksters", edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<edm::InputTag>("genP", edm::InputTag("genParticles"));
  // desc.add<edm::InputTag>("hitMap", edm::InputTag("hgcalRecHitMapProducer"));
  desc.add<edm::InputTag>("caloparticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("simCluster", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<double>("min_angle", 0.9);
  desc.add<bool>("do2Ddot", false);
  descriptions.add("tracksterLinkingDebugger", desc);
}
void TrackstersLinkingDebugger::beginRun(edm::Run const&, edm::EventSetup const& es) {
  const CaloGeometry& geom = es.getData(caloGeometry_token_);
  rhtools_.setGeometry(geom);
}
void TrackstersLinkingDebugger::beginJob() {
  // std::cout << "Minimum cosine " << min_angle_ << std::endl;
  printf("Minimum cosine %f Use 2D %d \n", min_angle_, do2D);
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "tree");
  tree_->Branch("cp_pt", &cp_pt);
  tree_->Branch("cp_E", &cp_E);
  tree_->Branch("simT_pt", &simT_pt);
  tree_->Branch("simT_E", &simT_E);
  tree_->Branch("candidate_pt", &candidate_pt);
  tree_->Branch("candidate_E", &candidate_E);
  tree_->Branch("fullMatch", &fullMatch);
  tree_->Branch("gen_pt", &gen_pt);
  tree_->Branch("gen_E", &gen_E);
  // tree_->Branch("simTrackstersEnergy", &simTrackstersEnergy);
  // tree_->Branch("simTrackstersRawEnergy", &simTrackstersRawEnergy);
}

void TrackstersLinkingDebugger::clear() {
  cp_pt->clear();
  cp_E->clear();
  simT_pt->clear();
  simT_E->clear();
  candidate_pt->clear();
  candidate_E->clear();
  gen_pt->clear();
  gen_E->clear();
  fullMatch->clear();
  // simTrackstersRawEnergy->clear();
  // genEnergySim->clear();
}

bool TrackstersLinkingDebugger::isCP(Trackster& t,
                                     CaloParticle& cp,
                                     edm::ProductID cp_id,
                                     edm::Handle<std::vector<CaloParticle>> cp_handle) {
  auto productID = t.seedID();
  auto index = t.seedIndex();
  bool result = false;
  if (productID == cp_id) {
    result = true;
    cp = (*cp_handle)[index];
  }
  return result;
}

void TrackstersLinkingDebugger::tracksterInfo(const Trackster& t, const edm::Event& event) const {
  auto e_over_h = (t.raw_em_pt() / ((t.raw_pt() - t.raw_em_pt()) != 0. ? (t.raw_pt() - t.raw_em_pt()) : 1.));
  std::cout << "\nTrackster raw_pt: " << t.raw_pt() << " raw_em_pt: " << t.raw_em_pt() << " eoh: " << e_over_h
            << " barycenter: " << t.barycenter() << " eta,phi (baricenter): " << t.barycenter().eta() << ", "
            << t.barycenter().phi() << " eta,phi (eigen): " << t.eigenvectors(0).eta() << ", "
            << t.eigenvectors(0).phi() << " pt(eigen): " << std::sqrt(t.eigenvectors(0).Unit().perp2()) * t.raw_energy()
            << " seedID: " << t.seedID() << " seedIndex: " << t.seedIndex() << " size: " << t.vertices().size()
            << " average usage: "
            << (std::accumulate(std::begin(t.vertex_multiplicity()), std::end(t.vertex_multiplicity()), 0.) /
                (float)t.vertex_multiplicity().size())
            << " raw_energy: " << t.raw_energy() << " regressed energy: " << t.regressed_energy()
            << " Vertex Multiplicity: " << (float)t.vertex_multiplicity().size() << " probs(ga/e/mu/np/cp/nh/am/unk): ";

  for (auto const& p : t.id_probabilities()) {
    std::cout << std::fixed << p << " ";
  }
  std::cout << " sigmas: ";
  for (auto const& s : t.sigmas()) {
    std::cout << s << " ";
  }
  std::cout << std::endl;
}

std::pair<std::vector<Trackster>, std::vector<int>> TrackstersLinkingDebugger::findCompatibleTrackster(
    std::vector<Trackster> tracksters,
    const double deltaR,
    const double min_cos,
    int tShift,
    int Nsize,
    std::vector<int>& mask) {
  std::vector<Trackster> compatibleTrackster;
  std::vector<int> compatibleTracksterId;
  std::pair<std::vector<Trackster>, std::vector<int>> results;
  auto t1 = tracksters[tShift];
  auto direction1 = t1.eigenvectors()[0].Unit();
  auto barycenter1 = t1.barycenter();
  auto eta1 = barycenter1.eta();
  auto phi1 = barycenter1.phi();

  auto howMany = 0.;
  for (int t_id2 = tShift + 1; t_id2 != Nsize; t_id2++) {
    if (tracksters[t_id2].raw_energy() == 0) {
      mask[t_id2] = 1;
      continue;
    }
    if (mask[t_id2] == 0) {
      auto t2 = tracksters[t_id2];
      auto direction2 = t2.eigenvectors()[0].Unit();
      auto barycenter2 = t2.barycenter();
      auto eta2 = barycenter2.eta();
      auto phi2 = barycenter2.phi();
      bool InDR = reco::deltaR(eta1, phi1, eta2, phi2) < deltaR ? true : false;

      if (InDR) {
        auto dot = direction1.Dot(direction2);

        if (dot > min_cos) {
          compatibleTrackster.emplace_back(t2);
          compatibleTracksterId.push_back(t_id2);
          mask[t_id2] = 1;
        }
        
      }
    }
  }
  // if (!compatibleTrackster.empty()) {
  compatibleTrackster.emplace_back(t1);
  compatibleTracksterId.push_back(tShift);
  mask[tShift] = 1;
  // }
  results = std::make_pair(compatibleTrackster, compatibleTracksterId);

  return results;
}

const int TrackstersLinkingDebugger::associateTrack(std::vector<Trackster> compatibleTrackster,
                                                    std::vector<reco::Track> tracks) {
  // For now return only if one of the linked trackster is seeded by a track.
  auto trackIdx = -99;
  for (auto t : compatibleTrackster) {
    if (t.seedIndex() >= 0) {
      trackIdx = t.seedIndex();
    }
  }
  // std::cout << "TRKIDX: " << trackIdx << std::endl;
  return trackIdx;
}

bool TrackstersLinkingDebugger::isNuclearInteraction(
    std::vector<Trackster> tracksters, double dR, int t_id1, int trackstersAvailable, std::vector<int>& mask) {
  auto countTrackster = 0;
  bool results = false;
  std::vector<int> tmpId;
  for (int t_id2 = t_id1 + 1; t_id2 != (int)(trackstersAvailable); t_id2++) {
    auto t1b = tracksters[t_id1].barycenter();
    auto t2b = tracksters[t_id2].barycenter();
    auto eta1 = t1b.eta();
    auto eta2 = t2b.eta();
    auto phi1 = t1b.phi();
    auto phi2 = t2b.phi();
    auto DR = reco::deltaR(eta1, phi1, eta2, phi2);
    tmpId.push_back(t_id2);
    if (DR < dR) {
      countTrackster += 1;
    }
  }
  if (countTrackster > 1) {
    for (auto id : tmpId) {
      mask[id] = 1;
    }
    mask[t_id1] = 1;
    results = true;
  }
  return results;
}

math::XYZTLorentzVector TrackstersLinkingDebugger::assignP4(Trackster t, std::vector<reco::Track> tracks) {
  math::XYZTLorentzVector p4Result;
  if (t.ticlIteration() == Trackster::IterationIndex::EM) {
    math::XYZTLorentzVector p4(t.raw_energy() * t.barycenter().unit().x(),
                               t.raw_energy() * t.barycenter().unit().y(),
                               t.raw_energy() * t.barycenter().unit().z(),
                               t.raw_energy());
    p4Result = p4;
  } else if (t.ticlIteration() == Trackster::IterationIndex::TRKEM ||
             t.ticlIteration() == Trackster::IterationIndex::TRKHAD) {
    auto seedIndex = t.seedIndex();
    auto track = tracks[seedIndex];
    math::XYZTLorentzVector p4(t.raw_energy() * track.momentum().unit().x(),
                               t.raw_energy() * track.momentum().unit().y(),
                               t.raw_energy() * track.momentum().unit().z(),
                               t.raw_energy());
    p4Result = p4;
  } else if (t.ticlIteration() == Trackster::IterationIndex::HAD) {
    constexpr double mpion = 0.13957;
    constexpr double mpion2 = mpion * mpion;
    float momentum = std::sqrt(t.raw_energy() * t.raw_energy() - mpion2);
    math::XYZTLorentzVector p4(momentum * t.barycenter().unit().x(),
                               momentum * t.barycenter().unit().y(),
                               momentum * t.barycenter().unit().z(),
                               t.raw_energy());
    p4Result = p4;
  }
  return p4Result;
}

void TrackstersLinkingDebugger::analyze(const edm::Event& evt, const edm::EventSetup& es) {
  // std::cout << evt.id() << "\n\n";

  evt_count += 1;
  std::cout << "Event " << evt_count << "\n" << std::endl;

  clear();
  const auto& tracksters = evt.get(tracksters_token_);
  const auto& simtracksters = evt.get(simtracksters_token_);
  const auto& layerClusters = evt.get(layerClusters_token);
  const auto& tracks = evt.get(tracks_token_);
  const auto& caloparticles = evt.get(caloparticles_token_);
  auto trackstersAvailable = (int)tracksters.size();
  std::vector<int> mask(trackstersAvailable, 0);
  std::vector<int> simMask(simtracksters.size(), 0);
  std::vector<int> trkMask(tracks.size(), 0);
  std::vector<Candidate> candidates;
  edm::Handle<std::vector<CaloParticle>> caloparticles_handle = evt.getHandle(caloparticles_token_);
  edm::Handle<std::vector<SimCluster>> simclusters_handle = evt.getHandle(clusters_token_);

  const auto& genParticles = evt.get(gensToken_);
  double pos = 0.;
  double neg = 0.;
  
  // double neg_zero = 0;
  for (auto tt : tracksters) {
    if (tt.raw_energy() > 0) {
      if (tt.barycenter().z() < 0) {
        neg += 1;
      } else if (tt.barycenter().z() > 0) {
        pos += 1;
      }
    }
  }
  if(neg == 1){
    tot_real_alone += 1;
  }
  if(pos == 1){
    tot_real_alone += 1;
  }

  for (int t_id1 = 0; t_id1 != trackstersAvailable - 1; t_id1++) {  // loop over trackster
    Candidate cand;
    if(mask[t_id1] == 0){
      
      auto compatibleTrackster = findCompatibleTrackster(tracksters, 1, min_angle_, t_id1, trackstersAvailable, mask);
      for(int c_id = 0; c_id != (int)compatibleTrackster.first.size(); c_id++){
        cand.linkedTrackster.push_back(compatibleTrackster.first[c_id]);
        cand.linkedTracksterId.push_back(compatibleTrackster.second[c_id]);
      }    
    }
    candidates.push_back(cand);
  }

  for(int id = 0; id != (int)mask.size(); id++){
    if(mask[id] == 0){
      Candidate candNotMask;
      candNotMask.linkedTrackster.push_back(tracksters[id]);
      candNotMask.linkedTracksterId.push_back(id);
      candidates.push_back(candNotMask);
    }
  }
  // Candidate candLast;
  // if(mask[trackstersAvailable-1] == 0 && tracksters[trackstersAvailable - 1].raw_energy() > 0){
  //   candLast.linkedTrackster.push_back(tracksters[trackstersAvailable-1]);
  //   candLast.linkedTracksterId.push_back(trackstersAvailable-1);
  //   candidates.push_back(candLast);
  // }

// Associate Track to candidate
  for(int c_id = 0; c_id != (int)candidates.size(); c_id++){
    auto candidate = candidates[c_id];
    if(!(candidate.linkedTracksterId.empty())){
      for(int t_id = 0; t_id != (int)candidate.linkedTrackster.size(); t_id++){
        auto t = candidate.linkedTrackster[t_id];
        if(t.seedIndex() >= 0){
          candidate.matchedTracks = tracks[t.seedIndex()];
          candidate.matchedTracksId = t.seedIndex();
          trkMask[t.seedIndex()] = 1;
        }
      }
      // if(candidate.matchedTracksId == -99){
      //   printf("Track not found %d\n", c_id);
      // }
    }
  }

//Comparison with SimTrackster
bool tmpZero = true;
  for(int sim_id = 0; sim_id != (int)simtracksters.size(); sim_id++){
    auto simZ = simtracksters[sim_id].barycenter().z();
    auto simT = simtracksters[sim_id];
    reco::GenParticle genp;
    CaloParticle cp;
    bool flagCP =
        isCP(simT, cp, caloparticles_handle.id(), caloparticles_handle);  //only if calo (no nuclear interaction)
    if(flagCP == true){ // no nuclear interaction
      if(tmpZero == true){
        for(auto ttmp : tracksters){
          if(pos == 0 || neg == 0){
            if(ttmp.raw_energy() == 0){
              zero += 1;
              tmpZero = false;
              printf("ZERO\n ");
            }
          }
        }
        
      }
    if(simZ < 0 && neg >= 2){
      real_double += 1;
    }
    if(simZ > 0 && pos >=2){
      real_double += 1;
    }
    tot_linkable_trackster += 1;
    printf("LINKABLE TRACKSTER \n");
    for(auto gP : genParticles){
      if(gP.eta()*cp.eta() > 0){
        genp = gP;
      }
    }
      for(int c_id = 0; c_id != (int)candidates.size(); c_id++){
        auto candidate = candidates[c_id];
        if(candidate.linkedTrackster.size() > 0){
          auto barycenterZ = candidate.linkedTrackster[0].barycenter().z();
          bool flagCompleteMatch = false;
          if(barycenterZ * simZ > 0){
            if(barycenterZ < 0){
              if(neg >= 1){
                printf("Percentage match %f ", (double)candidate.linkedTracksterId.size()/(double)neg);
                if((double)candidate.linkedTracksterId.size()/(double)neg == 1)
                  flagCompleteMatch = true;
                if(neg >= 2){
                  if(flagCompleteMatch == true){
                    tot_trackster_linked += 1;
                    printf("MATCHED ");
                  }
                  else{
                    tot_missed_link += 1;
                    printf("Missed ");
                  }
                }
              }
            }
            if(barycenterZ > 0){
              if(pos >= 1){
                printf("Percentage match %f ",(double)candidate.linkedTracksterId.size()/(double)pos);
                if((double)candidate.linkedTracksterId.size()/(double)pos == 1)
                  flagCompleteMatch = true;
                if(pos >= 2){
                  if(flagCompleteMatch == true){
                    printf("MATCHED ");
                    tot_trackster_linked += 1;
                  }
                  else{
                    tot_missed_link += 1;
                    printf("Missed ");
                  }
                }
              }
            }
             math::XYZTLorentzVector totP4(0,0,0,0);
             printf("Candidate %d Trackster ", c_id);
            for(int t_c_id = 0; t_c_id != (int)candidate.linkedTrackster.size(); t_c_id++){
              auto t_c = candidate.linkedTrackster[t_c_id];
              totP4 += assignP4(t_c, tracks);
              printf("%d ", candidate.linkedTracksterId[t_c_id]);
            }
            printf("CP %d pT %f Cp E %f TotP4: %f TotE: %f\n", sim_id,cp.pt(), cp.energy(), std::sqrt(totP4.perp2()), totP4.E());
            if(flagCompleteMatch){
              //fill variables;
              cp_pt->push_back(cp.pt());
              cp_E->push_back(cp.energy());
              simT_pt->push_back(simT.raw_pt());
              simT_E->push_back(simT.raw_energy());
              candidate_pt->push_back(std::sqrt(totP4.perp2()));
              candidate_E->push_back(totP4.E());
              gen_pt->push_back(genp.pt());
              gen_E->push_back(genp.energy());
              fullMatch->push_back(true);
            }
            else{
              cp_pt->push_back(0);
              cp_E->push_back(0);
              simT_pt->push_back(0);
              simT_E->push_back(0);
              candidate_pt->push_back(0);
              candidate_E->push_back(0);
              gen_pt->push_back(0);
              gen_E->push_back(0);
              fullMatch->push_back(false);

            }
          }
        }
      }
    }
  }

//Check candidates manually
  // for(auto c : candidates){
  //   printf("Candidates ");
  //   if(!(c.linkedTracksterId.empty())){
  //     for(auto id : c.linkedTracksterId){
  //       printf("%d ",id);
  //     }
  //   }
  //   printf("\n");
  // }
  tree_->Fill();
}

void TrackstersLinkingDebugger::endJob() {
  printf("Total CaloParticles %f\n", total_cp);
  printf("Linkable Trackster %f\n", tot_linkable_trackster);
  printf("Total Trackster Linked %f\n", tot_trackster_linked);
  printf("ZERO ENERGY %f \n", zero);
  printf("Total missed Link %f\n", tot_missed_link);  
  printf("Total Trackster Alone %f \n", tot_alone);
  printf("Total Trackster Zero %f \n", zeroT);
  // printf("Real Double %f \n", tot_linkable_trackster - (tot_real_alone + zero));
  printf("Real Double %f \n", real_double);
  printf("Real Alone %f \n", tot_real_alone + zero);
  printf("Tot Cand %f \n", tot_cand_found);
}
DEFINE_FWK_MODULE(TrackstersLinkingDebugger);
