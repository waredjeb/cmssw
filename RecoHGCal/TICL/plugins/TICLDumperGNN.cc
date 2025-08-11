#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <variant>
#include <algorithm>  // for std::minmax_element

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
#include "DataFormats/HGCalReco/interface/TICLGraph.h"
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

using CaloObjectVariant = std::variant<CaloParticle, SimCluster>;

namespace {
  template <typename TTMap>
  void calcReco2SimTracksterFit(const TTMap& recoToSimMap,
                                std::vector<int>& y,
                                std::vector<float>& sharedE,
                                std::vector<float>& score) {
    const std::size_t nReco = recoToSimMap.getMap().size();
    y.assign(nReco, -1);
    sharedE.assign(nReco, -1.f);
    score.assign(nReco, -1.f);

    for (std::size_t reco = 0; reco < nReco; ++reco) {
      std::size_t nPassed = 0;
      int candIdx = -1;
      float candSE = -1.f;
      float candSc = -1.f;

      for (const auto& sim : recoToSimMap[reco]) {
        if (sim.score() < 0.2f) {
          ++nPassed;
          candIdx = sim.index();
          candSE = sim.sharedEnergy();
          candSc = sim.score();
        }
      }
      if (nPassed == 1) {
        y[reco] = candIdx;
        sharedE[reco] = candSE;
        score[reco] = candSc;
      }
    }
  }

  inline void edgeLabelAndWeight(int src,
                                 int dst,
                                 const std::vector<int>& y,
                                 const std::vector<float>& delta,
                                 const std::vector<float>& sharedE,
                                 const std::vector<float>& rawE,
                                 int& label,
                                 float& weight) {
    const bool validEdge = (y[src] != -1) && (y[dst] != -1) && (y[src] == y[dst]);

    if (validEdge) {
      label = 1;

      float termSrc =
          (rawE[src] > std::numeric_limits<float>::epsilon()) ? (1.f - delta[src]) * sharedE[src] / rawE[src] : 0.f;
      float termDst =
          (rawE[dst] > std::numeric_limits<float>::epsilon()) ? (1.f - delta[dst]) * sharedE[dst] / rawE[dst] : 0.f;

      weight = (termSrc + termDst) / 2.;
    } else {
      label = 0;
      weight = 0.f;
    }
  }
}  // namespace

using TracksterToTracksterMap =
    ticl::AssociationMap<ticl::mapWithSharedEnergyAndScore, std::vector<ticl::Trackster>, std::vector<ticl::Trackster>>;
class TICLDumperGNN : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
public:
  float detector_size = (2 * (3 - 1.5) * (2 * 47));
  TICLDumperGNN(const edm::ParameterSet& params);
  ~TICLDumperGNN() override;

  void clearVariables();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override {};
  void endJob() override;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> sim_tracksters_token_;
  const edm::EDGetTokenT<TICLGraph> ticl_graph_token_;
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> layer_clusters_token_;
  edm::EDGetTokenT<TracksterToTracksterMap> associations_simToReco_token_;  ///< The tokens for each assocation
  edm::EDGetTokenT<TracksterToTracksterMap> associations_recoToSim_token_;
  const edm::EDGetTokenT<std::vector<SimCluster>> simclusters_token_;
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloparticles_token_;
  std::vector<float> node_raw_energy;
  std::vector<float> node_raw_em_energy;
  std::vector<float> node_barycenter_x;
  std::vector<float> node_barycenter_y;
  std::vector<float> node_barycenter_z;
  std::vector<float> node_barycenter_eta;
  std::vector<float> node_barycenter_phi;
  std::vector<float> node_eVector0_x;
  std::vector<float> node_eVector0_y;
  std::vector<float> node_eVector0_z;
  std::vector<float> node_EV1;
  std::vector<float> node_EV2;
  std::vector<float> node_EV3;
  std::vector<float> node_sigmaPCA1;
  std::vector<float> node_sigmaPCA2;
  std::vector<float> node_sigmaPCA3;
  std::vector<float> node_photon_prob;
  std::vector<float> node_electron_prob;
  std::vector<float> node_muon_prob;
  std::vector<float> node_neutral_pion_prob;
  std::vector<float> node_charged_hadron_prob;
  std::vector<float> node_neutral_hadron_prob;
  std::vector<float> node_z_min;
  std::vector<float> node_z_max;
  std::vector<float> node_time;
  std::vector<float> node_time_error;
  std::vector<float> node_LC_density;
  std::vector<float> node_trackster_density;
  std::vector<float> node_num_LCs;
  std::vector<float> node_num_hits;

  std::vector<float> simTrackster_raw_energy;
  std::vector<float> simTrackster_true_energy;
  std::vector<int> simTrackster_isPU;
  std::vector<int> simTrackster_pdgID;

  std::vector<int> node_match_idx;
  std::vector<float> node_match_sharedE;
  std::vector<float> node_match_score;

  std::vector<std::vector<float>> edge_raw_energy;
  std::vector<std::vector<float>> edge_barycenter_z;
  std::vector<std::vector<float>> edge_barycenter_xy;
  std::vector<std::vector<float>> edge_time;
  std::vector<std::vector<float>> edge_eigenvector0;
  std::vector<std::vector<unsigned int>> edgeIndex_out;
  std::vector<std::vector<unsigned int>> edgeIndex_in;

  std::vector<std::vector<int>> edge_label;
  std::vector<std::vector<float>> edge_weight;
  TTree* gnnTree;
};

TICLDumperGNN::TICLDumperGNN(edm::ParameterSet const& params)
    : tracksters_token_(consumes<std::vector<ticl::Trackster>>(params.getParameter<edm::InputTag>("tracksters"))),
      sim_tracksters_token_(
          consumes<std::vector<ticl::Trackster>>(params.getParameter<edm::InputTag>("simTracksters"))),
      ticl_graph_token_(consumes<TICLGraph>(params.getParameter<edm::InputTag>("ticlGraph"))),
      layer_clusters_token_(
          consumes<std::vector<reco::CaloCluster>>(params.getParameter<edm::InputTag>("layerClusters"))),
      associations_simToReco_token_(consumes<TracksterToTracksterMap>(params.getParameter<edm::InputTag>("simToReco"))),
      associations_recoToSim_token_(
          consumes<TracksterToTracksterMap>(params.getParameter<edm::InputTag>("recoToSim"))),
      simclusters_token_(consumes(params.getParameter<edm::InputTag>("simclusters"))),
      caloparticles_token_(consumes(params.getParameter<edm::InputTag>("caloparticles")))
{}

void TICLDumperGNN::clearVariables() {
  node_raw_energy.clear();
  node_raw_em_energy.clear();
  node_barycenter_x.clear();
  node_barycenter_y.clear();
  node_barycenter_z.clear();
  node_barycenter_eta.clear();
  node_barycenter_phi.clear();
  node_eVector0_x.clear();
  node_eVector0_y.clear();
  node_eVector0_z.clear();
  node_EV1.clear();
  node_EV2.clear();
  node_EV3.clear();
  node_sigmaPCA1.clear();
  node_sigmaPCA2.clear();
  node_sigmaPCA3.clear();
  node_photon_prob.clear();
  node_electron_prob.clear();
  node_muon_prob.clear();
  node_neutral_pion_prob.clear();
  node_charged_hadron_prob.clear();
  node_neutral_hadron_prob.clear();
  node_z_min.clear();
  node_z_max.clear();
  node_time.clear();
  node_time_error.clear();
  node_LC_density.clear();
  node_trackster_density.clear();
  node_num_LCs.clear();
  node_num_hits.clear();

  simTrackster_raw_energy.clear();
  simTrackster_true_energy.clear();
  simTrackster_isPU.clear();
  simTrackster_pdgID.clear();

  node_match_idx.clear();
  node_match_sharedE.clear();
  node_match_score.clear();
  edge_raw_energy.clear();
  edge_barycenter_z.clear();
  edge_barycenter_xy.clear();
  edge_time.clear();
  edge_eigenvector0.clear();
  edge_label.clear();
  edge_weight.clear();
  edgeIndex_out.clear();
  edgeIndex_in.clear();
}
TICLDumperGNN::~TICLDumperGNN() { clearVariables(); };
void TICLDumperGNN::endJob() {}

void TICLDumperGNN::beginJob() {
  edm::Service<TFileService> fs;
  gnnTree = fs->make<TTree>("GNNTraining", "GNNTraining");
  gnnTree->Branch("node_raw_energy", &node_raw_energy);
  gnnTree->Branch("node_raw_em_energy", &node_raw_em_energy);
  gnnTree->Branch("node_barycenter_x", &node_barycenter_x);
  gnnTree->Branch("node_barycenter_y", &node_barycenter_y);
  gnnTree->Branch("node_barycenter_z", &node_barycenter_z);
  gnnTree->Branch("node_barycenter_eta", &node_barycenter_eta);
  gnnTree->Branch("node_barycenter_phi", &node_barycenter_phi);
  gnnTree->Branch("node_eVector0_x", &node_eVector0_x);
  gnnTree->Branch("node_eVector0_y", &node_eVector0_y);
  gnnTree->Branch("node_eVector0_z", &node_eVector0_z);
  gnnTree->Branch("node_EV1", &node_EV1);
  gnnTree->Branch("node_EV2", &node_EV2);
  gnnTree->Branch("node_EV3", &node_EV3);
  gnnTree->Branch("node_sigmaPCA1", &node_sigmaPCA1);
  gnnTree->Branch("node_sigmaPCA2", &node_sigmaPCA2);
  gnnTree->Branch("node_sigmaPCA3", &node_sigmaPCA3);
  gnnTree->Branch("node_photon_prob", &node_photon_prob);
  gnnTree->Branch("node_electron_prob", &node_electron_prob);
  gnnTree->Branch("node_muon_prob", &node_muon_prob);
  gnnTree->Branch("node_neutral_pion_prob", &node_neutral_pion_prob);
  gnnTree->Branch("node_charged_hadron_prob", &node_charged_hadron_prob);
  gnnTree->Branch("node_neutral_hadron_prob", &node_neutral_hadron_prob);
  gnnTree->Branch("node_z_min", &node_z_min);
  gnnTree->Branch("node_z_max", &node_z_max);
  gnnTree->Branch("node_time", &node_time);
  gnnTree->Branch("node_time_error", &node_time_error);
  gnnTree->Branch("node_LC_density", &node_LC_density);
  gnnTree->Branch("node_trackster_density", &node_trackster_density);
  gnnTree->Branch("node_num_LCs", &node_num_LCs);
  gnnTree->Branch("node_num_hits", &node_num_hits);
  gnnTree->Branch("node_num_hits", &node_num_hits);
  gnnTree->Branch("node_match_idx", &node_match_idx);
  gnnTree->Branch("node_match_sharedE", &node_match_sharedE);
  gnnTree->Branch("node_match_score", &node_match_score);

  gnnTree->Branch("simTrackster_raw_energy", &simTrackster_raw_energy);
  gnnTree->Branch("simTrackster_true_energy", &simTrackster_true_energy);
  gnnTree->Branch("simTrackster_isPU", &simTrackster_isPU);
  gnnTree->Branch("simTrackster_pdgID", &simTrackster_pdgID);

  gnnTree->Branch("edge_raw_energy", &edge_raw_energy);
  gnnTree->Branch("edge_barycenter_z", &edge_barycenter_z);
  gnnTree->Branch("edge_barycenter_xy", &edge_barycenter_xy);
  gnnTree->Branch("edge_time", &edge_time);
  gnnTree->Branch("edge_eigenvector0", &edge_eigenvector0);
  gnnTree->Branch("edge_label", &edge_label);
  gnnTree->Branch("edge_weight", &edge_weight);
  gnnTree->Branch("edgeIndex_out", &edgeIndex_out);
  gnnTree->Branch("edgeIndex_in", &edgeIndex_in);
}
bool isFromPU(const ticl::Trackster& simTrackster,
              edm::Handle<std::vector<CaloParticle>>& caloparticles_h,
              const std::vector<CaloParticle>& caloparticles,
              const std::vector<SimCluster>& simclusters) {
  CaloObjectVariant caloObj;

  if (simTrackster.seedID() == caloparticles_h.id()) {
    caloObj = caloparticles[simTrackster.seedIndex()];

  } else {
    caloObj = simclusters[simTrackster.seedIndex()];
  }

  auto const& simTrack = std::visit([](auto&& obj) { return obj.g4Tracks()[0]; }, caloObj);

  if ((simTrack.eventId().event() != 0 or simTrack.eventId().bunchCrossing() != 0)) {
    return true;

  }

  else {
    return false;
  }
}
void TICLDumperGNN::beginRun(edm::Run const&, edm::EventSetup const& es) {};

void TICLDumperGNN::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  clearVariables();
  auto const& ticlGraph = event.get(ticl_graph_token_);
  auto const& tracksters = event.get(tracksters_token_);
  auto const& simTracksters = event.get(sim_tracksters_token_);
  auto const& simToRecoMap = event.get(associations_simToReco_token_);
  auto const& recoToSimMap = event.get(associations_recoToSim_token_);
  auto const& layerClusters = event.get(layer_clusters_token_);
  auto const& caloparticles = event.get(caloparticles_token_);
  auto caloparticles_h = event.getHandle(caloparticles_token_);
  auto simclusters_h = event.getHandle(simclusters_token_);
  auto const simclusters = event.get(simclusters_token_);

  // debug stream usage in concurrently scheduled modules

  int numTrackster = tracksters.size();
  int numEdges = ticlGraph.getNumberOfEdges();
  std::array<int, 3> const sizes{{numTrackster, numEdges, numEdges}};

  float trackster_dens = numTrackster / detector_size;

  std::vector<int> y;
  std::vector<float> sharedE, score;

  calcReco2SimTracksterFit(recoToSimMap, y, sharedE, score);

  node_match_idx = y;
  node_match_sharedE = sharedE;
  node_match_score = score;

  edge_label.resize(numTrackster);
  edge_weight.resize(numTrackster);

  for (size_t i = 0; i < simTracksters.size(); i++) {
    auto const& simT = simTracksters[i];
    simTrackster_raw_energy.push_back(simT.raw_energy());
    simTrackster_true_energy.push_back(simT.regressed_energy());
    CaloObjectVariant caloObj;
    if (simT.seedID() == caloparticles_h.id()) {
      caloObj = caloparticles[simT.seedIndex()];
    } else {
      caloObj = simclusters[simT.seedIndex()];
    }
    simTrackster_pdgID.push_back(std::visit([](auto&& obj) { return obj.pdgId(); }, caloObj));
    bool isPU = isFromPU(simT, caloparticles_h, caloparticles, simclusters);
    simTrackster_isPU.push_back(isPU);
  }
  for (int i = 0; i < numTrackster; i++) {
    node_raw_energy.push_back(tracksters[i].raw_energy());
    node_raw_em_energy.push_back(tracksters[i].raw_em_energy());
    node_barycenter_x.push_back(tracksters[i].barycenter().x());
    node_barycenter_y.push_back(tracksters[i].barycenter().y());
    node_barycenter_z.push_back(tracksters[i].barycenter().z());
    node_barycenter_eta.push_back(tracksters[i].barycenter().eta());
    node_barycenter_phi.push_back(tracksters[i].barycenter().phi());
    node_eVector0_x.push_back(tracksters[i].eigenvectors(0).x());
    node_eVector0_y.push_back(tracksters[i].eigenvectors(0).y());
    node_eVector0_z.push_back(tracksters[i].eigenvectors(0).z());
    node_EV1.push_back(tracksters[i].eigenvalues()[0]);
    node_EV2.push_back(tracksters[i].eigenvalues()[1]);
    node_EV3.push_back(tracksters[i].eigenvalues()[2]);
    node_sigmaPCA1.push_back(tracksters[i].sigmasPCA()[0]);
    node_sigmaPCA2.push_back(tracksters[i].sigmasPCA()[1]);
    node_sigmaPCA3.push_back(tracksters[i].sigmasPCA()[2]);
    node_photon_prob.push_back(tracksters[i].id_probability(ticl::Trackster::ParticleType::photon));
    node_electron_prob.push_back(tracksters[i].id_probability(ticl::Trackster::ParticleType::electron));
    node_muon_prob.push_back(tracksters[i].id_probability(ticl::Trackster::ParticleType::muon));
    node_neutral_pion_prob.push_back(tracksters[i].id_probability(ticl::Trackster::ParticleType::neutral_pion));
    node_charged_hadron_prob.push_back(tracksters[i].id_probability(ticl::Trackster::ParticleType::charged_hadron));
    node_neutral_hadron_prob.push_back(tracksters[i].id_probability(ticl::Trackster::ParticleType::neutral_hadron));
    node_num_LCs.push_back(tracksters[i].vertices().size());
    node_LC_density.push_back(tracksters[i].vertices().size() / detector_size);
    node_trackster_density.push_back(trackster_dens);
    node_time.push_back(tracksters[i].time());
    node_time_error.push_back(tracksters[i].timeError());

    int hits = 0;

    const auto& vertices = tracksters[i].vertices();

    auto minmax = std::minmax_element(vertices.begin(), vertices.end(), [&](const auto& a, const auto& b) {
      return layerClusters[a].z() < layerClusters[b].z();
    });

    node_z_min.push_back(layerClusters[*minmax.first].z());
    node_z_max.push_back(layerClusters[*minmax.second].z());

    for (const auto& vertex : vertices) {
      hits += layerClusters[vertex].size();
    }
    node_num_hits.push_back(hits);

    std::vector<unsigned int> outer = ticlGraph.getNode(i).getOuterNeighbours();
    edge_raw_energy.resize(numTrackster);
    edge_barycenter_z.resize(numTrackster);
    edge_barycenter_xy.resize(numTrackster);
    edge_time.resize(numTrackster);
    edge_eigenvector0.resize(numTrackster);
    edgeIndex_in.resize(numTrackster);
    edgeIndex_out.resize(numTrackster);
    edge_label.resize(numTrackster);
    edge_weight.resize(numTrackster);
  }

  for (int i = 0; i < numTrackster; i++) {
    std::vector<unsigned int> outer = ticlGraph.getNode(i).getOuterNeighbours();

    for (unsigned int node : outer) {
      edge_raw_energy[i].push_back(std::abs(tracksters[i].raw_energy() -
                                            tracksters[node].raw_energy()));  //todo: difference between raw energies
      edge_barycenter_z[i].push_back(std::abs(tracksters[i].barycenter().z() - tracksters[node].barycenter().z()));
      edge_time[i].push_back(std::abs(tracksters[i].time() - tracksters[node].time()));
      edge_barycenter_xy[i].push_back(std::hypot((tracksters[i].barycenter().x() - tracksters[node].barycenter().x()),
                                                 (tracksters[i].barycenter().y() - tracksters[node].barycenter().y())));

      edge_eigenvector0[i].push_back(
          std::acos(std::clamp((tracksters[i].eigenvectors(0).x() * tracksters[node].eigenvectors(0).x() +
                    tracksters[i].eigenvectors(0).y() * tracksters[node].eigenvectors(0).y() +
                    tracksters[i].eigenvectors(0).z() * tracksters[node].eigenvectors(0).z()),-1.f,1.f)));

      edgeIndex_in[i].push_back(i);
      edgeIndex_out[i].push_back(node);
      int edge_lab;
      float edge_w;

      edgeLabelAndWeight(i, node, y, score, sharedE, node_raw_energy, edge_lab, edge_w);
      edge_label[i].push_back(edge_lab);
      edge_weight[i].push_back(edge_w);
    }
  }
  gnnTree->Fill();
}

void TICLDumperGNN::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracksters", edm::InputTag("ticlTrackstersCLUE3DHigh"));
  desc.add<edm::InputTag>("ticlGraph", edm::InputTag("ticlGraph"));
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("simTracksters", edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>(
      "simToReco",
      edm::InputTag("allTrackstersToSimTrackstersAssociationsByLCs:ticlSimTrackstersToticlTrackstersCLUE3DHigh"));
  desc.add<edm::InputTag>(
      "recoToSim",
      edm::InputTag("allTrackstersToSimTrackstersAssociationsByLCs:ticlTrackstersCLUE3DHighToticlSimTracksters"));
  desc.add<edm::InputTag>("simclusters", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("caloparticles", edm::InputTag("mix", "MergedCaloTruth"));
  descriptions.add("ticlDumperGNN", desc);
}

DEFINE_FWK_MODULE(TICLDumperGNN);
