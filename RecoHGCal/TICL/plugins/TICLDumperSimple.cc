// Authors:  Philipp Zehetner, Wahid Redjeb, Aurora Perego, Felice Pantaleo

#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <variant>

#include <memory>  // unique_ptr
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"
#include "FWCore/ParameterSet/interface/allowedValues.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Provenance/interface/EventID.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Common.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"
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

#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using TracksterToTracksterMap =
    ticl::AssociationMap<ticl::mapWithSharedEnergyAndScore, std::vector<ticl::Trackster>, std::vector<ticl::Trackster>>;
// Helper class that dumps a single trackster collection (either tracksters or simTracksters)
class TracksterDumperHelperSimple {
public:
  enum class TracksterType {
    Trackster,       ///< Regular trackster (from RECO)
    SimTracksterCP,  ///< SimTrackster from CaloParticle
    SimTracksterSC   ///< SimTrackster from SimCluster
  };

  static TracksterType tracksterTypeFromString(std::string str) {
    if (str == "Trackster")
      return TracksterType::Trackster;
    if (str == "SimTracksterCP")
      return TracksterType::SimTracksterCP;
    if (str == "SimTracksterSC")
      return TracksterType::SimTracksterSC;
    throw std::runtime_error("TICLDumperSimple : TracksterDumperHelperSimple : Invalid trackster type " + str);
  }

  /** tracksterType : dtermines additional information that will be saved (calo truth information, track information) */
  TracksterDumperHelperSimple(TracksterType tracksterType = TracksterType::Trackster) : tracksterType_(tracksterType) {}

  /**
   * To be called once after tree creation. eventId_ should be a pointer to the EventID.
   * *Do not copy/move or resize vector holding object after calling this function*
  */
  void initTree(TTree* trackster_tree_, edm::EventID* eventId_) {
    trackster_tree_->Branch("raw_energy", &trackster_raw_energy);
    trackster_tree_->Branch("barycenter_x", &trackster_barycenter_x);
    trackster_tree_->Branch("barycenter_y", &trackster_barycenter_y);
    trackster_tree_->Branch("barycenter_z", &trackster_barycenter_z);
    trackster_tree_->Branch("barycenter_eta", &trackster_barycenter_eta);
    trackster_tree_->Branch("barycenter_phi", &trackster_barycenter_phi);
  }

  void clearVariables() {
    trackster_raw_energy.clear();
    trackster_barycenter_x.clear();
    trackster_barycenter_y.clear();
    trackster_barycenter_z.clear();
    trackster_barycenter_eta.clear();
    trackster_barycenter_phi.clear();
  }

  void fillFromEvent(std::vector<ticl::Trackster> const& tracksters,
                     edm::Handle<std::vector<SimCluster>> simClusters_h,
                     edm::Handle<std::vector<CaloParticle>> caloparticles_h) {
    nTracksters = tracksters.size();
    for (auto trackster_iterator = tracksters.begin(); trackster_iterator != tracksters.end(); ++trackster_iterator) {
      //per-trackster analysis
      trackster_raw_energy.push_back(trackster_iterator->raw_energy());
      trackster_barycenter_x.push_back(trackster_iterator->barycenter().x());
      trackster_barycenter_y.push_back(trackster_iterator->barycenter().y());
      trackster_barycenter_z.push_back(trackster_iterator->barycenter().z());
      trackster_barycenter_eta.push_back(trackster_iterator->barycenter().eta());
      trackster_barycenter_phi.push_back(trackster_iterator->barycenter().phi());
  }
  }


private:
  TracksterType tracksterType_;

  unsigned int nTracksters;
  unsigned int nClusters;
  std::vector<float> trackster_raw_energy;
  std::vector<float> trackster_barycenter_x;
  std::vector<float> trackster_barycenter_y;
  std::vector<float> trackster_barycenter_z;
  std::vector<float> trackster_barycenter_eta;
  std::vector<float> trackster_barycenter_phi;
};

// Helper class to dump a TracksterToSimTrackster association map (dumps recoToSim and simToReco at the same time)
class TracksterToSimTracksterAssociationHelperSimple {
public:
  /**
   * To be called once after tree creation. Output branches will be named prefix_recoToSim/simToReco_suffix_score/sharedE/
   * branchPrefix : for example tsCLUE3D. branchSuffix : usually one of SC or CP.
   * *Do not copy/move or resize vector holding object after calling this function*
  */
  void initTree(TTree* tree, std::string branchPrefix, std::string branchSuffix) {
    tree->Branch((branchPrefix + "_recoToSim_" + branchSuffix).c_str(), &recoToSim);
    tree->Branch((branchPrefix + "_recoToSim_" + branchSuffix + "_score").c_str(), &recoToSim_score);
    tree->Branch((branchPrefix + "_recoToSim_" + branchSuffix + "_sharedE").c_str(), &recoToSim_sharedE);
    tree->Branch((branchPrefix + "_simToReco_" + branchSuffix).c_str(), &simToReco);
    tree->Branch((branchPrefix + "_simToReco_" + branchSuffix + "_score").c_str(), &simToReco_score);
    tree->Branch((branchPrefix + "_simToReco_" + branchSuffix + "_sharedE").c_str(), &simToReco_sharedE);
  }

  void clearVariables() {
    recoToSim.clear();
    recoToSim_score.clear();
    recoToSim_sharedE.clear();
    simToReco.clear();
    simToReco_score.clear();
    simToReco_sharedE.clear();
  }

  void fillFromEvent(TracksterToTracksterMap const& recoToSimMap, TracksterToTracksterMap const& simToRecoMap) {
    // Reco -> Sim
    const auto numberOfTracksters = recoToSimMap.getMap().size();
    recoToSim.resize(numberOfTracksters);
    recoToSim_score.resize(numberOfTracksters);
    recoToSim_sharedE.resize(numberOfTracksters);

    for (size_t i = 0; i < numberOfTracksters; ++i) {
      for (const auto& simTracksterElement : recoToSimMap[i]) {
        recoToSim[i].push_back(simTracksterElement.index());
        recoToSim_sharedE[i].push_back(simTracksterElement.sharedEnergy());
        recoToSim_score[i].push_back(simTracksterElement.score());
      }
    }

    // Sim -> Reco
    const auto numberOfSimTracksters = simToRecoMap.getMap().size();
    simToReco.resize(numberOfSimTracksters);
    simToReco_score.resize(numberOfSimTracksters);
    simToReco_sharedE.resize(numberOfSimTracksters);

    for (size_t i = 0; i < numberOfSimTracksters; ++i) {
      for (const auto& recoTracksterElement : simToRecoMap[i]) {
        simToReco[i].push_back(recoTracksterElement.index());
        simToReco_sharedE[i].push_back(recoTracksterElement.sharedEnergy());
        simToReco_score[i].push_back(recoTracksterElement.score());
      }
    }
  }

private:
  std::vector<std::vector<uint32_t>> recoToSim;
  std::vector<std::vector<float>> recoToSim_score;
  std::vector<std::vector<float>> recoToSim_sharedE;
  std::vector<std::vector<uint32_t>> simToReco;
  std::vector<std::vector<float>> simToReco_score;
  std::vector<std::vector<float>> simToReco_sharedE;
};

class TICLDumperSimple : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
public:
  explicit TICLDumperSimple(const edm::ParameterSet&);
  ~TICLDumperSimple() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  typedef ticl::Vector Vector;
  typedef std::vector<double> Vec;

private:
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override {};
  void endJob() override;

  // Define Tokens
  const std::vector<edm::ParameterSet>
      tracksters_parameterSets_;  ///< A parameter set for each trackster collection to dump (giving tree name, etc)
  std::vector<edm::EDGetTokenT<std::vector<ticl::Trackster>>>
      tracksters_token_;  ///< a token for each trackster collection to dump
  std::vector<TracksterDumperHelperSimple>
      tracksters_dumperHelpers_;         ///< the dumper helpers for each trackster collection to dump
  std::vector<TTree*> tracksters_trees;  ///< TTree for each trackster collection to dump

  const std::vector<edm::ParameterSet>
      associations_parameterSets_;  ///< A parameter set for each associator collection to dump (with treeName, etc)
  std::vector<edm::EDGetTokenT<TracksterToTracksterMap>>
      associations_simToReco_token_;  ///< The tokens for each assocation
  std::vector<edm::EDGetTokenT<TracksterToTracksterMap>> associations_recoToSim_token_;
  std::vector<TracksterToSimTracksterAssociationHelperSimple>
      associations_dumperHelpers_;  ///< the dumper helpers for each association map to dump

  TTree* associations_tree_;

  const edm::EDGetTokenT<std::vector<SimCluster>> simclusters_token_;
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloparticles_token_;

  bool saveLCs_;
  bool saveSuperclustering_;
  bool saveSuperclusteringDNNScore_;
  bool saveRecoSuperclusters_;
  bool saveTICLCandidate_;
  bool saveSimTICLCandidate_;
  bool saveTracks_;

  // Output tree
  TTree* tree_;

  void clearVariables();

  // Variables for branches
  edm::EventID eventId_;
  unsigned int nclusters_;

  TTree* cluster_tree_;
  TTree* candidate_tree_;
  TTree* superclustering_tree_;
  TTree* tracks_tree_;
  TTree* simTICLCandidate_tree;
};

void TICLDumperSimple::clearVariables() {
  // event info
  nclusters_ = 0;

  for (TracksterDumperHelperSimple& tsDumper : tracksters_dumperHelpers_) {
    tsDumper.clearVariables();
  }

  for (auto& helper : associations_dumperHelpers_) {
    helper.clearVariables();
  }
};

TICLDumperSimple::TICLDumperSimple(const edm::ParameterSet& ps)
    : tracksters_parameterSets_(ps.getParameter<std::vector<edm::ParameterSet>>("tracksterCollections")),
      tracksters_token_(),
      associations_parameterSets_(ps.getParameter<std::vector<edm::ParameterSet>>("associators")),
      // The DumperHelpers should not be moved after construction (needed by TTree branch pointers), so construct them all here
      associations_dumperHelpers_(associations_parameterSets_.size()),
      simclusters_token_(consumes(ps.getParameter<edm::InputTag>("simclusters"))),
      caloparticles_token_(consumes(ps.getParameter<edm::InputTag>("caloparticles"))),
      saveLCs_(ps.getParameter<bool>("saveLCs")),
      saveSuperclustering_(ps.getParameter<bool>("saveSuperclustering")),
      //saveSuperclusteringDNNScore_(ps.getParameter<bool>("saveSuperclusteringDNNScore")),
      saveRecoSuperclusters_(ps.getParameter<bool>("saveRecoSuperclusters")),
      saveTICLCandidate_(ps.getParameter<bool>("saveSimTICLCandidate")),
      saveSimTICLCandidate_(ps.getParameter<bool>("saveSimTICLCandidate")),
      saveTracks_(ps.getParameter<bool>("saveTracks")) {
  for (edm::ParameterSet const& tracksterPset : tracksters_parameterSets_) {
    tracksters_token_.push_back(
        consumes<std::vector<ticl::Trackster>>(tracksterPset.getParameter<edm::InputTag>("inputTag")));
    tracksters_dumperHelpers_.emplace_back(
        TracksterDumperHelperSimple::tracksterTypeFromString(tracksterPset.getParameter<std::string>("tracksterType")));
  }

  for (edm::ParameterSet const& associationPset : associations_parameterSets_) {
    associations_recoToSim_token_.push_back(consumes<TracksterToTracksterMap>(
        edm::InputTag(associationPset.getParameter<edm::InputTag>("associatorRecoToSimInputTag"))));
    associations_simToReco_token_.push_back(consumes<TracksterToTracksterMap>(
        edm::InputTag(associationPset.getParameter<edm::InputTag>("associatorSimToRecoInputTag"))));
  }
};

TICLDumperSimple::~TICLDumperSimple() { clearVariables(); };

void TICLDumperSimple::beginRun(edm::Run const&, edm::EventSetup const& es) {
}

// Define tree and branches
void TICLDumperSimple::beginJob() {
  edm::Service<TFileService> fs;

  // Trackster trees
  for (unsigned int i = 0; i < tracksters_parameterSets_.size(); i++) {
    edm::ParameterSet const& tracksterPset = tracksters_parameterSets_[i];
    TTree* tree =
        fs->make<TTree>(tracksterPset.getParameter<std::string>("treeName").c_str(),
                        ("Tracksters : " + tracksterPset.getParameter<std::string>("treeName") +
                         " (InputTag : " + tracksterPset.getParameter<edm::InputTag>("inputTag").encode() + ")")
                            .c_str());
    tracksters_trees.push_back(tree);
    tracksters_dumperHelpers_[i].initTree(tree, &eventId_);
  }
  if (saveLCs_) {
    cluster_tree_ = fs->make<TTree>("clusters", "TICL tracksters");
  }
  if (saveTICLCandidate_) {
    candidate_tree_ = fs->make<TTree>("candidates", "TICL candidates");
  }
  if (saveSuperclustering_ || saveRecoSuperclusters_) {
    superclustering_tree_ = fs->make<TTree>("superclustering", "Superclustering in HGCAL CE-E");
  }
  if (saveRecoSuperclusters_) {
  }

  if (associations_parameterSets_.size() > 0) {
    associations_tree_ = fs->make<TTree>("associations", "Associations");
    associations_tree_->Branch("event", &eventId_);
  }
  for (unsigned int i = 0; i < associations_parameterSets_.size(); i++) {
    associations_dumperHelpers_[i].initTree(associations_tree_,
                                            associations_parameterSets_[i].getParameter<std::string>("branchName"),
                                            associations_parameterSets_[i].getParameter<std::string>("suffix"));
  }

  if (saveTracks_) {
    tracks_tree_ = fs->make<TTree>("tracks", "Tracks");
  }

  if (saveSimTICLCandidate_) {
    simTICLCandidate_tree = fs->make<TTree>("simTICLCandidate", "Sim TICL Candidate");
  }
}

void TICLDumperSimple::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  eventId_ = event.id();
  clearVariables();

  edm::Handle<std::vector<ticl::Trackster>> tracksters_in_candidate_handle;

  edm::Handle<std::vector<CaloParticle>> caloparticles_h;
  event.getByToken(caloparticles_token_, caloparticles_h);

  auto simclusters_h = event.getHandle(simclusters_token_);


  // Save all the trackster collections
  for (unsigned int i = 0; i < tracksters_dumperHelpers_.size(); i++) {
    edm::Handle<std::vector<ticl::Trackster>> tracksters_handle;
    std::vector<ticl::Trackster> const& tracksters = event.get<std::vector<ticl::Trackster>>(tracksters_token_[i]);
    tracksters_dumperHelpers_[i].fillFromEvent(
        tracksters, simclusters_h, caloparticles_h);
    tracksters_trees[i]->Fill();
  }

  // trackster to simTrackster associations
  for (unsigned int i = 0; i < associations_dumperHelpers_.size(); i++) {
    associations_dumperHelpers_[i].fillFromEvent(event.get(associations_recoToSim_token_[i]),
                                                 event.get(associations_simToReco_token_[i]));
  }
  if (associations_dumperHelpers_.size() > 0)
    associations_tree_->Fill();

  if (saveLCs_)
    cluster_tree_->Fill();
  if (saveTICLCandidate_)
    candidate_tree_->Fill();
  if (saveSuperclustering_ || saveRecoSuperclusters_)
    superclustering_tree_->Fill();
  if (saveTracks_)
    tracks_tree_->Fill();
  if (saveSimTICLCandidate_)
    simTICLCandidate_tree->Fill();
}

void TICLDumperSimple::endJob() {}

void TICLDumperSimple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  // Settings for dumping trackster collections
  edm::ParameterSetDescription tracksterDescValidator;
  tracksterDescValidator.add<std::string>("treeName")
      ->setComment("Name of the output tree for the trackster collection");
  tracksterDescValidator.add<edm::InputTag>("inputTag")->setComment("Input tag for the trackster collection to write");
  tracksterDescValidator.ifValue(
      edm::ParameterDescription<std::string>(
          "tracksterType",
          "Trackster",
          true,
          edm::Comment("Type of trackster. Trackster=regular trackster (from RECO). SimTracksterCP=Simtrackster "
                       "from CaloParticle. SimTracksterSC=Simtrackster from SimCluster")),
      edm::allowedValues<std::string>("Trackster", "SimTracksterCP", "SimTracksterSC"));
  desc.addVPSet("tracksterCollections", tracksterDescValidator)->setComment("Trackster collections to dump");

  // Settings for dumping trackster associators (recoToSim & simToReco)
  edm::ParameterSetDescription associatorDescValidator;
  associatorDescValidator.add<std::string>("branchName")->setComment("Name of the output branches in the tree");
  associatorDescValidator.add<std::string>("suffix")->setComment("Should be CP or SC (for the output branch name)");
  associatorDescValidator.add<edm::InputTag>("associatorRecoToSimInputTag")
      ->setComment("Input tag for the RecoToSim associator to dump");
  associatorDescValidator.add<edm::InputTag>("associatorSimToRecoInputTag")
      ->setComment("Input tag for the SimToReco associator to dump");
  desc.addVPSet("associators", associatorDescValidator)->setComment("Tracksters to SimTracksters associators to dump");

  desc.add<edm::InputTag>("simclusters", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("caloparticles", edm::InputTag("mix", "MergedCaloTruth"));

  desc.add<bool>("saveLCs", true);
  desc.add<bool>("saveTICLCandidate", true);
  desc.add<bool>("saveSimTICLCandidate", true);
  desc.add<bool>("saveTracks", true);
  desc.add<bool>("saveSuperclustering", true);
  desc.add<bool>("saveRecoSuperclusters", true)
      ->setComment("Save superclustering Egamma collections (as reco::SuperCluster)");
  descriptions.add("ticlDumperSimple", desc);
}

DEFINE_FWK_MODULE(TICLDumperSimple);
