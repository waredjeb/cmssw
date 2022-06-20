// -*- C++ -*-
//
// Package:    NTupler
// Class:      ClusterNtupler
//
/**\class RecHitTupler RecHitTupler.cc NTupler/plugins/ClusterNtupler.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Philipp Zehetner
// Based on example written by Thorben Quast
//
//


#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <sstream>

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

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "SimCalorimetry/HGCalAssociatorProducers/interface/AssociatorTools.h"
#include "SimDataFormats/Associations/interface/TracksterToSimTracksterHitLCAssociator.h"
#include "RecoHGCal/TICL/interface/commons.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

class Ntupler : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
    public:
  explicit Ntupler(const edm::ParameterSet&);
  ~Ntupler() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  typedef math::XYZVector Vector;
  typedef std::vector<double> Vec;
private:
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override{};
  void endJob() override;


  // some options
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_token_;
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> layer_clusters_token_;
  const edm::EDGetTokenT<TICLGraph> ticl_graph_token_;
  const edm::EDGetTokenT<std::vector<TICLCandidate>> ticl_candidates_token_;
  const edm::EDGetTokenT<std::vector<reco::Track>> tracks_token_;
  const edm::EDGetTokenT<std::vector<bool>> tracks_mask_token_;
  const edm::EDGetTokenT<edm::ValueMap<float>> tracks_time_token_;
  const edm::EDGetTokenT<edm::ValueMap<float>> tracks_time_quality_token_;
  const edm::EDGetTokenT<edm::ValueMap<float>> tracks_time_err_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_x_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_y_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_z_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_eta_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_phi_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_px_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_py_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_pz_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_merged_token_;
  const edm::EDGetTokenT<std::vector<float>> layerClustersLocalDensity_token_;
  const edm::EDGetTokenT<std::vector<float>> layerClustersRadius_token_;
  const edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> clustersTime_token_;
  const edm::EDGetTokenT<std::vector<int>> tracksterSeeds_token_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometry_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> simTracksters_SC_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> simTracksters_CP_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimSC_token_;
  const edm::EDGetTokenT<hgcal::SimToRecoCollectionSimTracksters> tsSimToRecoSC_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimCP_token_;
  const edm::EDGetTokenT<hgcal::SimToRecoCollectionSimTracksters> tsSimToRecoCP_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters> MergeRecoToSimSC_token_;
  const edm::EDGetTokenT<hgcal::SimToRecoCollectionSimTracksters> MergeSimToRecoSC_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters> MergeRecoToSimCP_token_;
  const edm::EDGetTokenT<hgcal::SimToRecoCollectionSimTracksters> MergeSimToRecoCP_token_;
  hgcal::RecHitTools rhtools_;
  

  // Output tree
  TTree* tree_;

  void clearVariables();

  unsigned int event_index;

  // Variables for branches
  unsigned int ev_event_;
  unsigned int ntracksters_;
  unsigned int stsSC_ntracksters_;
  unsigned int stsCP_ntracksters_;
  unsigned int nclusters_;
  unsigned int stsSC_nclusters_;
  unsigned int stsCP_nclusters_;
  size_t nsimTrackstersSC;
  size_t nsimTrackstersCP;

  std::vector<float_t> trackster_time;
  std::vector<float_t> trackster_timeError;
  std::vector<float_t> trackster_regressed_energy;
  std::vector<float_t> trackster_raw_energy;
  std::vector<float_t> trackster_raw_em_energy;
  std::vector<float_t> trackster_raw_pt;
  std::vector<float_t> trackster_raw_em_pt;
  std::vector<float_t> trackster_barycenter_x;
  std::vector<float_t> trackster_barycenter_y;
  std::vector<float_t> trackster_barycenter_z;
  std::vector<float_t> trackster_EV1;
  std::vector<float_t> trackster_EV2;
  std::vector<float_t> trackster_EV3;
  std::vector<float_t> trackster_eVector0_x;
  std::vector<float_t> trackster_eVector0_y;
  std::vector<float_t> trackster_eVector0_z;
  std::vector<float_t> trackster_sigmaPCA1;
  std::vector<float_t> trackster_sigmaPCA2;
  std::vector<float_t> trackster_sigmaPCA3;
  std::vector<float_t> trackster_barycenter_eta;
  std::vector<float_t> trackster_barycenter_phi;
  std::vector<std::vector<float_t>> trackster_id_probabilities;
  std::vector<std::vector<uint16_t> > trackster_vertices_indexes;
  std::vector<std::vector<float_t> > trackster_vertices_x;
  std::vector<std::vector<float_t> > trackster_vertices_y;
  std::vector<std::vector<float_t> > trackster_vertices_z;
  std::vector<std::vector<float_t> > trackster_vertices_energy;
  std::vector<std::vector<float_t> > trackster_vertices_correctedEnergy;
  std::vector<std::vector<float_t> > trackster_vertices_correctedEnergyUncertainty;
  std::vector<std::vector<float_t> > trackster_vertices_multiplicity;
  
  std::vector<float_t> stsSC_trackster_time;
  std::vector<float_t> stsSC_trackster_timeError;
  std::vector<float_t> stsSC_trackster_regressed_energy;
  std::vector<float_t> stsSC_trackster_raw_energy;
  std::vector<float_t> stsSC_trackster_raw_em_energy;
  std::vector<float_t> stsSC_trackster_raw_pt;
  std::vector<float_t> stsSC_trackster_raw_em_pt;
  std::vector<float_t> stsSC_trackster_barycenter_x;
  std::vector<float_t> stsSC_trackster_barycenter_y;
  std::vector<float_t> stsSC_trackster_barycenter_z;
  std::vector<float_t> stsSC_trackster_EV1;
  std::vector<float_t> stsSC_trackster_EV2;
  std::vector<float_t> stsSC_trackster_EV3;
  std::vector<float_t> stsSC_trackster_eVector0_x;
  std::vector<float_t> stsSC_trackster_eVector0_y;
  std::vector<float_t> stsSC_trackster_eVector0_z;
  std::vector<float_t> stsSC_trackster_sigmaPCA1;
  std::vector<float_t> stsSC_trackster_sigmaPCA2;
  std::vector<float_t> stsSC_trackster_sigmaPCA3;
  std::vector<float_t> stsSC_trackster_barycenter_eta;
  std::vector<float_t> stsSC_trackster_barycenter_phi;
  std::vector<std::vector<float_t>> stsSC_trackster_id_probabilities;
  std::vector<std::vector<uint16_t> > stsSC_trackster_vertices_indexes;
  std::vector<std::vector<float_t> > stsSC_trackster_vertices_x;
  std::vector<std::vector<float_t> > stsSC_trackster_vertices_y;
  std::vector<std::vector<float_t> > stsSC_trackster_vertices_z;
  std::vector<std::vector<float_t> > stsSC_trackster_vertices_energy;
  std::vector<std::vector<float_t> > stsSC_trackster_vertices_correctedEnergy;
  std::vector<std::vector<float_t> > stsSC_trackster_vertices_correctedEnergyUncertainty;
  std::vector<std::vector<float_t> > stsSC_trackster_vertices_multiplicity;
  std::vector<float_t> stsCP_trackster_time;
  std::vector<float_t> stsCP_trackster_timeError;
  std::vector<float_t> stsCP_trackster_regressed_energy;
  std::vector<float_t> stsCP_trackster_raw_energy;
  std::vector<float_t> stsCP_trackster_raw_em_energy;
  std::vector<float_t> stsCP_trackster_raw_pt;
  std::vector<float_t> stsCP_trackster_raw_em_pt;
  std::vector<float_t> stsCP_trackster_barycenter_x;
  std::vector<float_t> stsCP_trackster_barycenter_y;
  std::vector<float_t> stsCP_trackster_barycenter_z;
  std::vector<float_t> stsCP_trackster_EV1;
  std::vector<float_t> stsCP_trackster_EV2;
  std::vector<float_t> stsCP_trackster_EV3;
  std::vector<float_t> stsCP_trackster_eVector0_x;
  std::vector<float_t> stsCP_trackster_eVector0_y;
  std::vector<float_t> stsCP_trackster_eVector0_z;
  std::vector<float_t> stsCP_trackster_sigmaPCA1;
  std::vector<float_t> stsCP_trackster_sigmaPCA2;
  std::vector<float_t> stsCP_trackster_sigmaPCA3;
  std::vector<float_t> stsCP_trackster_barycenter_eta;
  std::vector<float_t> stsCP_trackster_barycenter_phi;
  std::vector<std::vector<float_t>> stsCP_trackster_id_probabilities;
  std::vector<std::vector<uint16_t> > stsCP_trackster_vertices_indexes;
  std::vector<std::vector<float_t> > stsCP_trackster_vertices_x;
  std::vector<std::vector<float_t> > stsCP_trackster_vertices_y;
  std::vector<std::vector<float_t> > stsCP_trackster_vertices_z;
  std::vector<std::vector<float_t> > stsCP_trackster_vertices_energy;
  std::vector<std::vector<float_t> > stsCP_trackster_vertices_correctedEnergy;
  std::vector<std::vector<float_t> > stsCP_trackster_vertices_correctedEnergyUncertainty;
  std::vector<std::vector<float_t> > stsCP_trackster_vertices_multiplicity;

  
  std::vector<int> layer_cluster_seed;
  std::vector<int> stsSC_layer_cluster_seed;
  std::vector<int> stsCP_layer_cluster_seed;
  
  // from TICLGraph
  std::vector<std::vector<uint16_t>> node_linked_inners;
  std::vector<std::vector<uint16_t>> node_linked_outers;
  std::vector<bool> isRootTrackster;

  // from TICLCandidate, product of linking
  size_t nCandidates;
  std::vector<int> candidate_charge;
  std::vector<int> candidate_pdgId;
  std::vector<float_t> candidate_energy;
  std::vector<double> candidate_px;
  std::vector<double> candidate_py;
  std::vector<double> candidate_pz;
  std::vector<float_t> candidate_time;
  std::vector<float_t> candidate_time_err;
  std::vector<std::vector<float_t>> candidate_id_probabilities;
  std::vector<std::vector<uint16_t>> tracksters_in_candidate;
  std::vector<uint16_t> track_in_candidate;

  // merged tracksters
  size_t nTrackstersMerged;
  std::vector<float_t> tracksters_merged_barycenter_x;
  std::vector<float_t> tracksters_merged_barycenter_y;
  std::vector<float_t> tracksters_merged_barycenter_z;
  std::vector<float_t> tracksters_merged_EV1;
  std::vector<float_t> tracksters_merged_EV2;
  std::vector<float_t> tracksters_merged_EV3;
  std::vector<float_t> tracksters_merged_eVector0_x;
  std::vector<float_t> tracksters_merged_eVector0_y;
  std::vector<float_t> tracksters_merged_eVector0_z;
  std::vector<float_t> tracksters_merged_sigmaPCA1;
  std::vector<float_t> tracksters_merged_sigmaPCA2;
  std::vector<float_t> tracksters_merged_sigmaPCA3;
  std::vector<std::vector<float_t>> tracksters_merged_id_probabilities;

  // associations
  std::vector<std::vector<uint16_t>> trackstersCLUE3D_recoToSim_SC;
  std::vector<std::vector<float_t>> trackstersCLUE3D_recoToSim_SC_score;
  std::vector<std::vector<uint16_t>> trackstersCLUE3D_simToReco_SC;
  std::vector<std::vector<float_t>> trackstersCLUE3D_simToReco_SC_score;
  std::vector<std::vector<float_t>> trackstersCLUE3D_simToReco_SC_sharedE;

  std::vector<std::vector<uint16_t>> trackstersCLUE3D_recoToSim_CP;
  std::vector<std::vector<float_t>> trackstersCLUE3D_recoToSim_CP_score;
  std::vector<std::vector<uint16_t>> trackstersCLUE3D_simToReco_CP;
  std::vector<std::vector<float_t>> trackstersCLUE3D_simToReco_CP_score;
  std::vector<std::vector<float_t>> trackstersCLUE3D_simToReco_CP_sharedE;

  std::vector<std::vector<uint16_t>> MergeTracksters_recoToSim_SC;
  std::vector<std::vector<float_t>> MergeTracksters_recoToSim_SC_score;
  std::vector<std::vector<uint16_t>> MergeTracksters_simToReco_SC;
  std::vector<std::vector<float_t>> MergeTracksters_simToReco_SC_score;
  std::vector<std::vector<float_t>> MergeTracksters_simToReco_SC_sharedE;

  std::vector<std::vector<uint16_t>> MergeTracksters_recoToSim_CP;
  std::vector<std::vector<float_t>> MergeTracksters_recoToSim_CP_score;
  std::vector<std::vector<uint16_t>> MergeTracksters_simToReco_CP;
  std::vector<std::vector<float_t>> MergeTracksters_simToReco_CP_score;
  std::vector<std::vector<float_t>> MergeTracksters_simToReco_CP_sharedE;

  std::vector<uint32_t> cluster_seedID;
  std::vector<float_t> cluster_energy;
  std::vector<float_t> cluster_correctedEnergy;
  std::vector<float_t> cluster_correctedEnergyUncertainty;
  std::vector<float_t> cluster_position_x;
  std::vector<float_t> cluster_position_y;
  std::vector<float_t> cluster_position_z;
  std::vector<float_t> cluster_position_eta;
  std::vector<float_t> cluster_position_phi;
  std::vector<int> cluster_type;
  std::vector<float> cluster_time;
  std::vector<float> cluster_timeErr;
  std::vector<float> cluster_ld;
  std::vector<float> cluster_radius;
  std::vector<uint32_t> cluster_number_of_hits;


  std::vector<int>track_ev;
  std::vector<unsigned int>track_id;
  std::vector<float>track_hgcal_x;
  std::vector<float>track_hgcal_y;
  std::vector<float>track_hgcal_z;
  std::vector<float>track_hgcal_px;
  std::vector<float>track_hgcal_py;
  std::vector<float>track_hgcal_pz;
  std::vector<float>track_hgcal_eta;
  std::vector<float>track_hgcal_phi;
  std::vector<float>track_pt;
  std::vector<int>track_charge;
  std::vector<double>track_time;
  std::vector<float>track_time_quality;
  std::vector<float>track_time_err;
  std::vector<int>track_nhits;

  TTree* trackster_tree_;
  TTree* cluster_tree_;
  TTree* graph_tree_;
  TTree* candidate_tree_;
  TTree* tracksters_merged_tree_;
  TTree* associations_tree_;
  TTree* simtrackstersSC_tree_;
  TTree* simtrackstersCP_tree_;
  TTree* tracks_tree_;;
};

void Ntupler::clearVariables() {
  // event info
  ev_event_ = 0;
  ntracksters_ = 0;
  nclusters_ = 0;

  trackster_time.clear();
  trackster_timeError.clear();
  trackster_regressed_energy.clear();
  trackster_raw_energy.clear();
  trackster_raw_em_energy.clear();
  trackster_raw_pt.clear();
  trackster_raw_em_pt.clear();
  trackster_barycenter_x.clear();
  trackster_barycenter_y.clear();
  trackster_barycenter_z.clear();
  trackster_EV1.clear();
  trackster_EV2.clear();
  trackster_EV3.clear();
  trackster_eVector0_x.clear();
  trackster_eVector0_y.clear();
  trackster_eVector0_z.clear();
  trackster_sigmaPCA1.clear();
  trackster_sigmaPCA2.clear();
  trackster_sigmaPCA3.clear();
  trackster_barycenter_eta.clear();
  trackster_barycenter_phi.clear();
  trackster_id_probabilities.clear();
  trackster_vertices_indexes.clear();
  trackster_vertices_x.clear();
  trackster_vertices_y.clear();
  trackster_vertices_z.clear();
  trackster_vertices_energy.clear();
  trackster_vertices_correctedEnergy.clear();
  trackster_vertices_correctedEnergyUncertainty.clear();
  layer_cluster_seed.clear();

  stsSC_trackster_time.clear();
  stsSC_trackster_timeError.clear();
  stsSC_trackster_regressed_energy.clear();
  stsSC_trackster_raw_energy.clear();
  stsSC_trackster_raw_em_energy.clear();
  stsSC_trackster_raw_pt.clear();
  stsSC_trackster_raw_em_pt.clear();
  stsSC_trackster_barycenter_x.clear();
  stsSC_trackster_barycenter_y.clear();
  stsSC_trackster_barycenter_z.clear();
  stsSC_trackster_EV1.clear();
  stsSC_trackster_EV2.clear();
  stsSC_trackster_EV3.clear();
  stsSC_trackster_eVector0_x.clear();
  stsSC_trackster_eVector0_y.clear();
  stsSC_trackster_eVector0_z.clear();
  stsSC_trackster_sigmaPCA1.clear();
  stsSC_trackster_sigmaPCA2.clear();
  stsSC_trackster_sigmaPCA3.clear();
  stsSC_trackster_barycenter_eta.clear();
  stsSC_trackster_barycenter_phi.clear();
  stsSC_trackster_id_probabilities.clear();
  stsSC_trackster_vertices_indexes.clear();
  stsSC_trackster_vertices_x.clear();
  stsSC_trackster_vertices_y.clear();
  stsSC_trackster_vertices_z.clear();
  stsSC_trackster_vertices_energy.clear();
  stsSC_trackster_vertices_correctedEnergy.clear();
  stsSC_trackster_vertices_correctedEnergyUncertainty.clear();
  stsSC_layer_cluster_seed.clear();

  stsCP_trackster_time.clear();
  stsCP_trackster_timeError.clear();
  stsCP_trackster_regressed_energy.clear();
  stsCP_trackster_raw_energy.clear();
  stsCP_trackster_raw_em_energy.clear();
  stsCP_trackster_raw_pt.clear();
  stsCP_trackster_raw_em_pt.clear();
  stsCP_trackster_barycenter_x.clear();
  stsCP_trackster_barycenter_y.clear();
  stsCP_trackster_barycenter_z.clear();
  stsCP_trackster_EV1.clear();
  stsCP_trackster_EV2.clear();
  stsCP_trackster_EV3.clear();
  stsCP_trackster_eVector0_x.clear();
  stsCP_trackster_eVector0_y.clear();
  stsCP_trackster_eVector0_z.clear();
  stsCP_trackster_sigmaPCA1.clear();
  stsCP_trackster_sigmaPCA2.clear();
  stsCP_trackster_sigmaPCA3.clear();
  stsCP_trackster_barycenter_eta.clear();
  stsCP_trackster_barycenter_phi.clear();
  stsCP_trackster_id_probabilities.clear();
  stsCP_trackster_vertices_indexes.clear();
  stsCP_trackster_vertices_x.clear();
  stsCP_trackster_vertices_y.clear();
  stsCP_trackster_vertices_z.clear();
  stsCP_trackster_vertices_energy.clear();
  stsCP_trackster_vertices_correctedEnergy.clear();
  stsCP_trackster_vertices_correctedEnergyUncertainty.clear();
  stsCP_layer_cluster_seed.clear();

  node_linked_inners.clear();
  node_linked_outers.clear();
  isRootTrackster.clear();
  
  nCandidates = 0;
  candidate_charge.clear();
  candidate_pdgId.clear();
  candidate_energy.clear();
  candidate_px.clear();
  candidate_py.clear();
  candidate_pz.clear();
  candidate_time.clear();
  candidate_time_err.clear();
  candidate_id_probabilities.clear();
  tracksters_in_candidate.clear();
  track_in_candidate.clear();

  nTrackstersMerged = 0;
  tracksters_merged_barycenter_x.clear();
  tracksters_merged_barycenter_y.clear();
  tracksters_merged_barycenter_z.clear();
  tracksters_merged_EV1.clear();
  tracksters_merged_EV2.clear();
  tracksters_merged_EV3.clear();
  tracksters_merged_eVector0_x.clear();
  tracksters_merged_eVector0_y.clear();
  tracksters_merged_eVector0_z.clear();
  tracksters_merged_sigmaPCA1.clear();
  tracksters_merged_sigmaPCA2.clear();
  tracksters_merged_sigmaPCA3.clear();
  tracksters_merged_id_probabilities.clear();

  trackstersCLUE3D_recoToSim_SC.clear();
  trackstersCLUE3D_recoToSim_SC_score.clear();
  trackstersCLUE3D_simToReco_SC.clear();
  trackstersCLUE3D_simToReco_SC_score.clear();
  trackstersCLUE3D_simToReco_SC_sharedE.clear();

  trackstersCLUE3D_recoToSim_CP.clear();
  trackstersCLUE3D_recoToSim_CP_score.clear();
  trackstersCLUE3D_simToReco_CP.clear();
  trackstersCLUE3D_simToReco_CP_score.clear();
  trackstersCLUE3D_simToReco_CP_sharedE.clear();

  MergeTracksters_recoToSim_SC.clear();
  MergeTracksters_recoToSim_SC_score.clear();
  MergeTracksters_simToReco_SC.clear();
  MergeTracksters_simToReco_SC_score.clear();
  MergeTracksters_simToReco_SC_sharedE.clear();

  MergeTracksters_recoToSim_CP.clear();
  MergeTracksters_recoToSim_CP_score.clear();
  MergeTracksters_simToReco_CP.clear();
  MergeTracksters_simToReco_CP_score.clear();
  MergeTracksters_simToReco_CP_sharedE.clear();

  nsimTrackstersSC = 0;

  cluster_seedID.clear();
  cluster_energy.clear();
  cluster_correctedEnergy.clear();
  cluster_correctedEnergyUncertainty.clear();
  cluster_position_x.clear();
  cluster_position_y.clear();
  cluster_position_z.clear();
  cluster_position_eta.clear();
  cluster_position_phi.clear();
  cluster_type.clear();
  cluster_time.clear();
  cluster_timeErr.clear();
  cluster_ld.clear();
  cluster_radius.clear();
  cluster_number_of_hits.clear();

  track_ev.clear();
  track_id.clear();
  track_hgcal_x.clear();
  track_hgcal_y.clear();
  track_hgcal_z.clear();
  track_hgcal_eta.clear();
  track_hgcal_phi.clear();
  track_hgcal_px.clear();
  track_hgcal_py.clear();
  track_hgcal_pz.clear();
  track_pt.clear();
  track_charge.clear();
  track_time.clear();
  track_time_quality.clear();
  track_time_err.clear();
  track_nhits.clear();
};


Ntupler::Ntupler(const edm::ParameterSet& ps)
    : tracksters_token_(consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("trackstersclue3d"))),
      layer_clusters_token_(
          consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layerClusters"))),
      ticl_graph_token_(consumes<TICLGraph>(ps.getParameter<edm::InputTag>("ticlgraph"))),
      ticl_candidates_token_(consumes<std::vector<TICLCandidate>>(ps.getParameter<edm::InputTag>("ticlcandidates"))),
      tracks_token_(consumes<std::vector<reco::Track>>(ps.getParameter<edm::InputTag>("tracks"))),
      tracks_mask_token_(consumes<std::vector<bool>>(ps.getParameter<edm::InputTag>("masked_tracks"))),
      tracks_time_token_(consumes<edm::ValueMap<float>>(ps.getParameter<edm::InputTag>("tracksTime"))),
      tracks_time_quality_token_(consumes<edm::ValueMap<float>>(ps.getParameter<edm::InputTag>("tracksTimeQual"))),
      tracks_time_err_token_(consumes<edm::ValueMap<float>>(ps.getParameter<edm::InputTag>("tracksTimeErr"))),
      hgcaltracks_x_token_(consumes<std::vector<double>>(ps.getParameter<edm::InputTag>("hgcaltracks_x"))),
      hgcaltracks_y_token_(consumes<std::vector<double>>(ps.getParameter<edm::InputTag>("hgcaltracks_y"))),
      hgcaltracks_z_token_(consumes<std::vector<double>>(ps.getParameter<edm::InputTag>("hgcaltracks_z"))),
      hgcaltracks_eta_token_(consumes<std::vector<double>>(ps.getParameter<edm::InputTag>("hgcaltracks_eta"))),
      hgcaltracks_phi_token_(consumes<std::vector<double>>(ps.getParameter<edm::InputTag>("hgcaltracks_phi"))),
      hgcaltracks_px_token_(consumes<std::vector<double>>(ps.getParameter<edm::InputTag>("hgcaltracks_px"))),
      hgcaltracks_py_token_(consumes<std::vector<double>>(ps.getParameter<edm::InputTag>("hgcaltracks_py"))),
      hgcaltracks_pz_token_(consumes<std::vector<double>>(ps.getParameter<edm::InputTag>("hgcaltracks_pz"))),
      tracksters_merged_token_(consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("trackstersmerged"))),
      layerClustersLocalDensity_token_(consumes<std::vector<float>>(ps.getParameter<edm::InputTag>("layerClustersLocalDensity"))),
      layerClustersRadius_token_(consumes<std::vector<float>>(ps.getParameter<edm::InputTag>("layerClustersRadius"))),
      clustersTime_token_(
          consumes<edm::ValueMap<std::pair<float, float>>>(ps.getParameter<edm::InputTag>("layer_clustersTime"))),
      tracksterSeeds_token_(consumes<std::vector<int>>(ps.getParameter<edm::InputTag>("tracksterSeeds"))),
      caloGeometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>()),
      simTracksters_SC_token_(consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("simtrackstersSC"))),
      simTracksters_CP_token_(consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("simtrackstersCP"))),
      tsRecoToSimSC_token_(consumes<hgcal::RecoToSimCollectionSimTracksters>(ps.getParameter<edm::InputTag>("recoToSimAssociatorSC"))),
      tsSimToRecoSC_token_(consumes<hgcal::SimToRecoCollectionSimTracksters>(ps.getParameter<edm::InputTag>("simToRecoAssociatorSC"))),
      tsRecoToSimCP_token_(consumes<hgcal::RecoToSimCollectionSimTracksters>(ps.getParameter<edm::InputTag>("recoToSimAssociatorCP"))),
      tsSimToRecoCP_token_(consumes<hgcal::SimToRecoCollectionSimTracksters>(ps.getParameter<edm::InputTag>("simToRecoAssociatorCP"))),
      MergeRecoToSimSC_token_(consumes<hgcal::RecoToSimCollectionSimTracksters>(ps.getParameter<edm::InputTag>("MergerecoToSimAssociatorSC"))),
      MergeSimToRecoSC_token_(consumes<hgcal::SimToRecoCollectionSimTracksters>(ps.getParameter<edm::InputTag>("MergesimToRecoAssociatorSC"))),
      MergeRecoToSimCP_token_(consumes<hgcal::RecoToSimCollectionSimTracksters>(ps.getParameter<edm::InputTag>("MergerecoToSimAssociatorCP"))),
      MergeSimToRecoCP_token_(consumes<hgcal::SimToRecoCollectionSimTracksters>(ps.getParameter<edm::InputTag>("MergesimToRecoAssociatorCP")))      {
      };

Ntupler::~Ntupler() { clearVariables(); };

void Ntupler::beginRun(edm::Run const&, edm::EventSetup const& es) {
  const CaloGeometry& geom = es.getData(caloGeometry_token_);
  rhtools_.setGeometry(geom);
// Define tree and branches  
}


void Ntupler::beginJob() {
  edm::Service<TFileService> fs;
  trackster_tree_ = fs->make<TTree>("tracksters", "TICL tracksters");
  cluster_tree_ = fs->make<TTree>("clusters", "TICL tracksters");
  graph_tree_ = fs->make<TTree>("graph", "TICL graph");
  candidate_tree_ = fs->make<TTree>("candidates", "TICL candidates");
  tracksters_merged_tree_ = fs->make<TTree>("trackstersMerged", "TICL tracksters merged");
  associations_tree_ = fs->make<TTree>("associations", "Associations");
  simtrackstersSC_tree_ = fs->make<TTree>("simtrackstersSC", "TICL simTracksters SC");
  simtrackstersCP_tree_ = fs->make<TTree>("simtrackstersCP", "TICL simTracksters CP");
  tracks_tree_ = fs->make<TTree>("tracks", "Tracks");

  trackster_tree_->Branch("event", &ev_event_);
  trackster_tree_->Branch("NClusters", &nclusters_);
  trackster_tree_->Branch("NTracksters", &ntracksters_);
  trackster_tree_->Branch("time", &trackster_time);
  trackster_tree_->Branch("timeError", &trackster_timeError);
  trackster_tree_->Branch("regressed_energy", &trackster_regressed_energy);
  trackster_tree_->Branch("raw_energy", &trackster_raw_energy);
  trackster_tree_->Branch("raw_em_energy", &trackster_raw_em_energy);
  trackster_tree_->Branch("raw_pt", &trackster_raw_pt);
  trackster_tree_->Branch("raw_em_pt", &trackster_raw_em_pt);
  trackster_tree_->Branch("barycenter_x", &trackster_barycenter_x);
  trackster_tree_->Branch("barycenter_y", &trackster_barycenter_y);
  trackster_tree_->Branch("barycenter_z", &trackster_barycenter_z);
  trackster_tree_->Branch("trackster_barycenter_eta", &trackster_barycenter_eta);
  trackster_tree_->Branch("trackster_barycenter_phi", &trackster_barycenter_phi);
  trackster_tree_->Branch("EV1", &trackster_EV1);
  trackster_tree_->Branch("EV2", &trackster_EV2);
  trackster_tree_->Branch("EV3", &trackster_EV3);
  trackster_tree_->Branch("eVector0_x", &trackster_eVector0_x);
  trackster_tree_->Branch("eVector0_y", &trackster_eVector0_y);
  trackster_tree_->Branch("eVector0_z", &trackster_eVector0_z);
  trackster_tree_->Branch("sigmaPCA1", &trackster_sigmaPCA1);
  trackster_tree_->Branch("sigmaPCA2", &trackster_sigmaPCA2);
  trackster_tree_->Branch("sigmaPCA3", &trackster_sigmaPCA3);
  trackster_tree_->Branch("id_probabilities", &trackster_id_probabilities);
  trackster_tree_->Branch("vertices_indexes", &trackster_vertices_indexes);
  trackster_tree_->Branch("vertices_x", &trackster_vertices_x);
  trackster_tree_->Branch("vertices_y", &trackster_vertices_y);
  trackster_tree_->Branch("vertices_z", &trackster_vertices_z);
  trackster_tree_->Branch("vertices_energy", &trackster_vertices_energy);
  trackster_tree_->Branch("vertices_correctedEnergy", &trackster_vertices_correctedEnergy);
  trackster_tree_->Branch("vertices_correctedEnergyUncertainty", &trackster_vertices_correctedEnergyUncertainty);
  trackster_tree_->Branch("vertices_multiplicity", &trackster_vertices_multiplicity); //NEW
  trackster_tree_->Branch("layer_cluster_seed", &layer_cluster_seed);

  simtrackstersSC_tree_->Branch("stsSC_event", &ev_event_);
  simtrackstersSC_tree_->Branch("stsSC_NClusters", &stsSC_nclusters_);
  simtrackstersSC_tree_->Branch("stsSC_NTracksters", &stsSC_ntracksters_);
  simtrackstersSC_tree_->Branch("stsSC_time", &stsSC_trackster_time);
  simtrackstersSC_tree_->Branch("stsSC_timeError", &stsSC_trackster_timeError);
  simtrackstersSC_tree_->Branch("stsSC_regressed_energy", &stsSC_trackster_regressed_energy);
  simtrackstersSC_tree_->Branch("stsSC_raw_energy", &stsSC_trackster_raw_energy);
  simtrackstersSC_tree_->Branch("stsSC_raw_em_energy", &stsSC_trackster_raw_em_energy);
  simtrackstersSC_tree_->Branch("stsSC_raw_pt", &stsSC_trackster_raw_pt);
  simtrackstersSC_tree_->Branch("stsSC_raw_em_pt", &stsSC_trackster_raw_em_pt);
  simtrackstersSC_tree_->Branch("stsSC_barycenter_x", &stsSC_trackster_barycenter_x);
  simtrackstersSC_tree_->Branch("stsSC_barycenter_y", &stsSC_trackster_barycenter_y);
  simtrackstersSC_tree_->Branch("stsSC_barycenter_z", &stsSC_trackster_barycenter_z);
  simtrackstersSC_tree_->Branch("stsSC_trackster_barycenter_eta", &stsSC_trackster_barycenter_eta);
  simtrackstersSC_tree_->Branch("stsSC_trackster_barycenter_phi", &stsSC_trackster_barycenter_phi);
  simtrackstersSC_tree_->Branch("stsSC_EV1", &stsSC_trackster_EV1);
  simtrackstersSC_tree_->Branch("stsSC_EV2", &stsSC_trackster_EV2);
  simtrackstersSC_tree_->Branch("stsSC_EV3", &stsSC_trackster_EV3);
  simtrackstersSC_tree_->Branch("stsSC_eVector0_x", &stsSC_trackster_eVector0_x);
  simtrackstersSC_tree_->Branch("stsSC_eVector0_y", &stsSC_trackster_eVector0_y);
  simtrackstersSC_tree_->Branch("stsSC_eVector0_z", &stsSC_trackster_eVector0_z);
  simtrackstersSC_tree_->Branch("stsSC_sigmaPCA1", &stsSC_trackster_sigmaPCA1);
  simtrackstersSC_tree_->Branch("stsSC_sigmaPCA2", &stsSC_trackster_sigmaPCA2);
  simtrackstersSC_tree_->Branch("stsSC_sigmaPCA3", &stsSC_trackster_sigmaPCA3);
  simtrackstersSC_tree_->Branch("stsSC_trackster_barycenter_eta", &stsSC_trackster_barycenter_eta);
  simtrackstersSC_tree_->Branch("stsSC_trackster_barycenter_phi", &stsSC_trackster_barycenter_phi);
  simtrackstersSC_tree_->Branch("stsSC_id_probabilities", &stsSC_trackster_id_probabilities);
  simtrackstersSC_tree_->Branch("stsSC_vertices_indexes", &stsSC_trackster_vertices_indexes);
  simtrackstersSC_tree_->Branch("stsSC_vertices_x", &stsSC_trackster_vertices_x);
  simtrackstersSC_tree_->Branch("stsSC_vertices_y", &stsSC_trackster_vertices_y);
  simtrackstersSC_tree_->Branch("stsSC_vertices_z", &stsSC_trackster_vertices_z);
  simtrackstersSC_tree_->Branch("stsSC_vertices_energy", &stsSC_trackster_vertices_energy);
  simtrackstersSC_tree_->Branch("stsSC_vertices_correctedEnergy", &stsSC_trackster_vertices_correctedEnergy);
  simtrackstersSC_tree_->Branch("stsSC_vertices_correctedEnergyUncertainty", &stsSC_trackster_vertices_correctedEnergyUncertainty);
  simtrackstersSC_tree_->Branch("stsSC_vertices_multiplicity", &stsSC_trackster_vertices_multiplicity); //NEW
  simtrackstersSC_tree_->Branch("NsimTrackstersSC", &nsimTrackstersSC);

  simtrackstersCP_tree_->Branch("stsCP_event", &ev_event_);
  simtrackstersCP_tree_->Branch("stsCP_NClusters", &stsCP_nclusters_);
  simtrackstersCP_tree_->Branch("stsCP_NTracksters", &stsCP_ntracksters_);
  simtrackstersCP_tree_->Branch("stsCP_time", &stsCP_trackster_time);
  simtrackstersCP_tree_->Branch("stsCP_timeError", &stsCP_trackster_timeError);
  simtrackstersCP_tree_->Branch("stsCP_regressed_energy", &stsCP_trackster_regressed_energy);
  simtrackstersCP_tree_->Branch("stsCP_raw_energy", &stsCP_trackster_raw_energy);
  simtrackstersCP_tree_->Branch("stsCP_raw_em_energy", &stsCP_trackster_raw_em_energy);
  simtrackstersCP_tree_->Branch("stsCP_raw_pt", &stsCP_trackster_raw_pt);
  simtrackstersCP_tree_->Branch("stsCP_raw_em_pt", &stsCP_trackster_raw_em_pt);
  simtrackstersCP_tree_->Branch("stsCP_barycenter_x", &stsCP_trackster_barycenter_x);
  simtrackstersCP_tree_->Branch("stsCP_barycenter_y", &stsCP_trackster_barycenter_y);
  simtrackstersCP_tree_->Branch("stsCP_barycenter_z", &stsCP_trackster_barycenter_z);
  simtrackstersCP_tree_->Branch("stsCP_trackster_barycenter_eta", &stsCP_trackster_barycenter_eta);
  simtrackstersCP_tree_->Branch("stsCP_trackster_barycenter_phi", &stsCP_trackster_barycenter_phi);
  simtrackstersCP_tree_->Branch("stsCP_EV1", &stsCP_trackster_EV1);
  simtrackstersCP_tree_->Branch("stsCP_EV2", &stsCP_trackster_EV2);
  simtrackstersCP_tree_->Branch("stsCP_EV3", &stsCP_trackster_EV3);
  simtrackstersCP_tree_->Branch("stsCP_eVector0_x", &stsCP_trackster_eVector0_x);
  simtrackstersCP_tree_->Branch("stsCP_eVector0_y", &stsCP_trackster_eVector0_y);
  simtrackstersCP_tree_->Branch("stsCP_eVector0_z", &stsCP_trackster_eVector0_z);
  simtrackstersCP_tree_->Branch("stsCP_sigmaPCA1", &stsCP_trackster_sigmaPCA1);
  simtrackstersCP_tree_->Branch("stsCP_sigmaPCA2", &stsCP_trackster_sigmaPCA2);
  simtrackstersCP_tree_->Branch("stsCP_sigmaPCA3", &stsCP_trackster_sigmaPCA3);
  simtrackstersCP_tree_->Branch("stsCP_id_probabilities", &stsCP_trackster_id_probabilities);
  simtrackstersCP_tree_->Branch("stsCP_vertices_indexes", &stsCP_trackster_vertices_indexes);
  simtrackstersCP_tree_->Branch("stsCP_vertices_x", &stsCP_trackster_vertices_x);
  simtrackstersCP_tree_->Branch("stsCP_vertices_y", &stsCP_trackster_vertices_y);
  simtrackstersCP_tree_->Branch("stsCP_vertices_z", &stsCP_trackster_vertices_z);
  simtrackstersCP_tree_->Branch("stsCP_vertices_energy", &stsCP_trackster_vertices_energy);
  simtrackstersCP_tree_->Branch("stsCP_vertices_correctedEnergy", &stsCP_trackster_vertices_correctedEnergy);
  simtrackstersCP_tree_->Branch("stsCP_vertices_correctedEnergyUncertainty", &stsCP_trackster_vertices_correctedEnergyUncertainty);
  simtrackstersCP_tree_->Branch("stsCP_vertices_multiplicity", &stsCP_trackster_vertices_multiplicity); //NEW
  simtrackstersSC_tree_->Branch("NsimTrackstersCP", &nsimTrackstersCP);

  graph_tree_->Branch("linked_inners", &node_linked_inners);
  graph_tree_->Branch("linked_outers", &node_linked_outers);
  graph_tree_->Branch("isRootTrackster", &isRootTrackster);

  candidate_tree_->Branch("NCandidates", &nCandidates);
  candidate_tree_->Branch("candidate_charge", &candidate_charge);
  candidate_tree_->Branch("candidate_pdgId", &candidate_pdgId);
  candidate_tree_->Branch("candidate_id_probabilities", &candidate_id_probabilities);
  candidate_tree_->Branch("candidate_time", &candidate_time);
  candidate_tree_->Branch("candidate_timeErr", &candidate_time_err);
  candidate_tree_->Branch("candidate_energy", &candidate_energy);
  candidate_tree_->Branch("candidate_px", &candidate_px);
  candidate_tree_->Branch("candidate_py", &candidate_py);
  candidate_tree_->Branch("candidate_pz", &candidate_pz);
  candidate_tree_->Branch("track_in_candidate", &track_in_candidate);
  candidate_tree_->Branch("tracksters_in_candidate", &tracksters_in_candidate);

  tracksters_merged_tree_->Branch("NTrackstersMerged", &nTrackstersMerged);
  tracksters_merged_tree_->Branch("barycenter_x", &tracksters_merged_barycenter_x);
  tracksters_merged_tree_->Branch("barycenter_y", &tracksters_merged_barycenter_y);
  tracksters_merged_tree_->Branch("barycenter_z", &tracksters_merged_barycenter_z);
  tracksters_merged_tree_->Branch("EV1", &tracksters_merged_EV1);
  tracksters_merged_tree_->Branch("EV2", &tracksters_merged_EV2);
  tracksters_merged_tree_->Branch("EV3", &tracksters_merged_EV3);
  tracksters_merged_tree_->Branch("eVector0_x", &tracksters_merged_eVector0_x);
  tracksters_merged_tree_->Branch("eVector0_y", &tracksters_merged_eVector0_y);
  tracksters_merged_tree_->Branch("eVector0_z", &tracksters_merged_eVector0_z);
  tracksters_merged_tree_->Branch("sigmaPCA1", &tracksters_merged_sigmaPCA1);
  tracksters_merged_tree_->Branch("sigmaPCA2", &tracksters_merged_sigmaPCA2);
  tracksters_merged_tree_->Branch("sigmaPCA3", &tracksters_merged_sigmaPCA3);
  tracksters_merged_tree_->Branch("id_probabilities", &tracksters_merged_id_probabilities);

  associations_tree_->Branch("tsCLUE3D_recoToSim_SC", &trackstersCLUE3D_recoToSim_SC);
  associations_tree_->Branch("tsCLUE3D_recoToSim_SC_score", &trackstersCLUE3D_recoToSim_SC_score);
  associations_tree_->Branch("tsCLUE3D_simToReco_SC", &trackstersCLUE3D_simToReco_SC);
  associations_tree_->Branch("tsCLUE3D_simToReco_SC_score", &trackstersCLUE3D_simToReco_SC_score);
  associations_tree_->Branch("tsCLUE3D_recoToSim_SC_sharedE", &trackstersCLUE3D_simToReco_SC_sharedE);

  associations_tree_->Branch("tsCLUE3D_recoToSim_CP", &trackstersCLUE3D_recoToSim_CP);
  associations_tree_->Branch("tsCLUE3D_recoToSim_CP_score", &trackstersCLUE3D_recoToSim_CP_score);
  associations_tree_->Branch("tsCLUE3D_simToReco_CP", &trackstersCLUE3D_simToReco_CP);
  associations_tree_->Branch("tsCLUE3D_simToReco_CP_score", &trackstersCLUE3D_simToReco_CP_score);
  associations_tree_->Branch("tsCLUE3D_recoToSim_CP_sharedE", &trackstersCLUE3D_simToReco_CP_sharedE);

  associations_tree_->Branch("Mergetstracksters_recoToSim_SC", &MergeTracksters_recoToSim_SC);
  associations_tree_->Branch("Mergetstracksters_recoToSim_SC_score", &MergeTracksters_recoToSim_SC_score);
  associations_tree_->Branch("Mergetstracksters_simToReco_SC", &MergeTracksters_simToReco_SC);
  associations_tree_->Branch("Mergetstracksters_simToReco_SC_score", &MergeTracksters_simToReco_SC_score);
  associations_tree_->Branch("Mergetstracksters_recoToSim_SC_sharedE", &MergeTracksters_simToReco_SC_sharedE);

  associations_tree_->Branch("Mergetracksters_recoToSim_CP", &MergeTracksters_recoToSim_CP);
  associations_tree_->Branch("Mergetracksters_recoToSim_CP_score", &MergeTracksters_recoToSim_CP_score);
  associations_tree_->Branch("Mergetracksters_simToReco_CP", &MergeTracksters_simToReco_CP);
  associations_tree_->Branch("Mergetracksters_simToReco_CP_score", &MergeTracksters_simToReco_CP_score);
  associations_tree_->Branch("Mergetracksters_recoToSim_CP_sharedE", &MergeTracksters_simToReco_CP_sharedE);
  
  cluster_tree_->Branch("seedID", &cluster_seedID);
  cluster_tree_->Branch("energy", &cluster_energy);
  cluster_tree_->Branch("correctedEnergy", &cluster_correctedEnergy);
  cluster_tree_->Branch("correctedEnergyUncertainty", &cluster_correctedEnergyUncertainty);
  cluster_tree_->Branch("position_x", &cluster_position_x);
  cluster_tree_->Branch("position_y", &cluster_position_y);
  cluster_tree_->Branch("position_z", &cluster_position_z);
  cluster_tree_->Branch("position_eta", &cluster_position_eta);
  cluster_tree_->Branch("position_phi", &cluster_position_phi);
  cluster_tree_->Branch("cluster_type", &cluster_type);
  cluster_tree_->Branch("cluster_time", &cluster_time);
  cluster_tree_->Branch("cluster_timeErr", &cluster_timeErr);
  cluster_tree_->Branch("cluster_local_density", &cluster_ld);
  cluster_tree_->Branch("cluster_radius", &cluster_radius);
  cluster_tree_->Branch("cluster_number_of_hits", &cluster_number_of_hits);

  tracks_tree_->Branch("track_ev",&track_ev);
  tracks_tree_->Branch("track_id",&track_id);
  tracks_tree_->Branch("track_hgcal_x",&track_hgcal_x);
  tracks_tree_->Branch("track_hgcal_y",&track_hgcal_y);
  tracks_tree_->Branch("track_hgcal_z",&track_hgcal_z);
  tracks_tree_->Branch("track_hgcal_py",&track_hgcal_eta);
  tracks_tree_->Branch("track_hgcal_pz",&track_hgcal_phi);
  tracks_tree_->Branch("track_hgcal_px",&track_hgcal_px);
  tracks_tree_->Branch("track_hgcal_py",&track_hgcal_py);
  tracks_tree_->Branch("track_hgcal_pz",&track_hgcal_pz);
  tracks_tree_->Branch("track_pt",&track_pt);
  tracks_tree_->Branch("track_charge",&track_charge);
  tracks_tree_->Branch("track_time",&track_time);
  tracks_tree_->Branch("track_time_quality",&track_time_quality);
  tracks_tree_->Branch("track_time_err",&track_time_err);
  tracks_tree_->Branch("track_nhits",&track_nhits); 

  event_index = 0;
}

void Ntupler::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  event_index++;
  clearVariables();

  //get all the tracksters
  edm::Handle<std::vector<ticl::Trackster>> tracksters_handle;
  event.getByToken(tracksters_token_, tracksters_handle);
  const auto& tracksters = *tracksters_handle;

  //get all the layer clusters
  edm::Handle<std::vector<reco::CaloCluster>> layer_clusters_h;
  event.getByToken(layer_clusters_token_, layer_clusters_h);
  const auto& clusters = *layer_clusters_h;

  edm::Handle<edm::ValueMap<std::pair<float, float>>> clustersTime_h;
  event.getByToken(clustersTime_token_, clustersTime_h);
  const auto &layerClustersTimes = *clustersTime_h;

  //TICL Graph
  edm::Handle<TICLGraph> ticl_graph_h;
  event.getByToken(ticl_graph_token_, ticl_graph_h);
  const auto& graph = *ticl_graph_h;

  //TICL Candidate
  edm::Handle<std::vector<TICLCandidate>> candidates_h;
  event.getByToken(ticl_candidates_token_, candidates_h);
  const auto& ticlcandidates = *candidates_h;

  //Track
  edm::Handle<std::vector<reco::Track>> tracks_h;
  event.getByToken(tracks_token_, tracks_h);
  const auto& tracks = *tracks_h;

  edm::Handle<std::vector<bool>> mask_tracks_h;
  event.getByToken(tracks_mask_token_, mask_tracks_h);
  const auto& mask_tracks = *mask_tracks_h;

  edm::Handle<std::vector<double>> hgcaltracks_x_h;
  event.getByToken(hgcaltracks_x_token_, hgcaltracks_x_h);
  const auto& tracks_propagated_x = *hgcaltracks_x_h;

  edm::Handle<std::vector<double>> hgcaltracks_y_h;
  event.getByToken(hgcaltracks_y_token_, hgcaltracks_y_h);
  const auto& tracks_propagated_y = *hgcaltracks_y_h;

  edm::Handle<std::vector<double>> hgcaltracks_z_h;
  event.getByToken(hgcaltracks_z_token_, hgcaltracks_z_h);
  const auto& tracks_propagated_z = *hgcaltracks_z_h;

  edm::Handle<std::vector<double>> hgcaltracks_eta_h;
  event.getByToken(hgcaltracks_eta_token_, hgcaltracks_eta_h);
  const auto& tracks_propagated_eta = *hgcaltracks_eta_h;

  edm::Handle<std::vector<double>> hgcaltracks_phi_h;
  event.getByToken(hgcaltracks_phi_token_, hgcaltracks_phi_h);
  const auto& tracks_propagated_phi = *hgcaltracks_phi_h;
  
  edm::Handle<std::vector<double>> hgcaltracks_px_h;
  event.getByToken(hgcaltracks_px_token_, hgcaltracks_px_h);
  const auto& tracks_propagated_px = *hgcaltracks_px_h;

  edm::Handle<std::vector<double>> hgcaltracks_py_h;
  event.getByToken(hgcaltracks_py_token_, hgcaltracks_py_h);
  const auto& tracks_propagated_py = *hgcaltracks_py_h;

  edm::Handle<std::vector<double>> hgcaltracks_pz_h;
  event.getByToken(hgcaltracks_pz_token_, hgcaltracks_pz_h);
  const auto& tracks_propagated_pz = *hgcaltracks_pz_h;

  edm::Handle<edm::ValueMap<float>> trackTime_h;
  event.getByToken(tracks_time_token_, trackTime_h);
  const auto &trackTime = *trackTime_h;

  edm::Handle<edm::ValueMap<float>> trackTimeErr_h;
  event.getByToken(tracks_time_err_token_, trackTimeErr_h);
  const auto &trackTimeErr = *trackTimeErr_h;

  edm::Handle<edm::ValueMap<float>> trackTimeQual_h;
  event .getByToken(tracks_time_quality_token_, trackTimeQual_h);
  const auto &trackTimeQual = *trackTimeQual_h;

  //Tracksters merged
  edm::Handle<std::vector<ticl::Trackster>> tracksters_merged_h;
  event.getByToken(tracksters_merged_token_, tracksters_merged_h);
  const auto& trackstersmerged = *tracksters_merged_h;

  // simTracksters from SC
  edm::Handle<std::vector<ticl::Trackster>> simTrackstersSC_h;
  event.getByToken(simTracksters_SC_token_, simTrackstersSC_h);
  const auto& simTrackstersSC = *simTrackstersSC_h;

// simTracksters from CP
  edm::Handle<std::vector<ticl::Trackster>> simTrackstersCP_h;
  event.getByToken(simTracksters_CP_token_, simTrackstersCP_h);
  const auto& simTrackstersCP = *simTrackstersCP_h;

  // trackster reco to sim SC
  edm::Handle<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimSC_h;
  event.getByToken(tsRecoToSimSC_token_, tsRecoToSimSC_h);
  auto const& tsRecoSimSCMap = *tsRecoToSimSC_h;

  // sim simTrackster SC to reco trackster
  edm::Handle<hgcal::SimToRecoCollectionSimTracksters> tsSimToRecoSC_h;
  event.getByToken(tsSimToRecoSC_token_, tsSimToRecoSC_h);
  auto const& tsSimToRecoSCMap = *tsSimToRecoSC_h;

  // trackster reco to sim CP
  edm::Handle<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimCP_h;
  event.getByToken(tsRecoToSimCP_token_, tsRecoToSimCP_h);
  auto const& tsRecoSimCPMap = *tsRecoToSimCP_h;

  // sim simTrackster CP to reco trackster
  edm::Handle<hgcal::SimToRecoCollectionSimTracksters> tsSimToRecoCP_h;
  event.getByToken(tsSimToRecoCP_token_, tsSimToRecoCP_h);
  auto const& tsSimToRecoCPMap = *tsSimToRecoCP_h;

//
  edm::Handle<hgcal::RecoToSimCollectionSimTracksters> mergetsRecoToSimSC_h;
  event.getByToken(MergeRecoToSimSC_token_, mergetsRecoToSimSC_h);
  auto const& MergetsRecoSimSCMap = *mergetsRecoToSimSC_h;

  // sim simTrackster SC to reco trackster
  edm::Handle<hgcal::SimToRecoCollectionSimTracksters> mergetsSimToRecoSC_h;
  event.getByToken(MergeSimToRecoSC_token_, mergetsSimToRecoSC_h);
  auto const& MergetsSimToRecoSCMap = *mergetsSimToRecoSC_h;

  // trackster reco to sim CP
  edm::Handle<hgcal::RecoToSimCollectionSimTracksters> mergetsRecoToSimCP_h;
  event.getByToken(MergeRecoToSimCP_token_, mergetsRecoToSimCP_h);
  auto const& MergetsRecoSimCPMap = *mergetsRecoToSimCP_h;

  // sim simTrackster CP to reco trackster
  edm::Handle<hgcal::SimToRecoCollectionSimTracksters> mergetsSimToRecoCP_h;
  event.getByToken(MergeSimToRecoCP_token_, mergetsSimToRecoCP_h);
  auto const& MergetsSimToRecoCPMap = *mergetsSimToRecoCP_h;


  edm::Handle<std::vector<float>> layer_cluster_local_density_h;
  event.getByToken(layerClustersLocalDensity_token_, layer_cluster_local_density_h);
  const auto& layer_cluster_density = *layer_cluster_local_density_h;

  edm::Handle<std::vector<float>> layer_cluster_radius_h;
  event.getByToken(layerClustersRadius_token_, layer_cluster_radius_h);
  const auto& layer_cluster_radius = *layer_cluster_radius_h;

  edm::Handle<std::vector<int>> trackster_cluster_seed_h;
  event.getByToken(tracksterSeeds_token_, trackster_cluster_seed_h);
  const auto& tracksterSeeds = *trackster_cluster_seed_h;



  ev_event_ = event_index;
  ntracksters_ = tracksters.size();
  nclusters_ = clusters.size();

  int t_id = 0;
  for (auto trackster_iterator = tracksters.begin(); trackster_iterator != tracksters.end(); ++trackster_iterator) {
      //per-trackster analysis
    trackster_time.push_back(trackster_iterator->time());
    trackster_timeError.push_back(trackster_iterator->timeError());
    trackster_regressed_energy.push_back(trackster_iterator->regressed_energy());
    trackster_raw_energy.push_back(trackster_iterator->raw_energy());
    trackster_raw_em_energy.push_back(trackster_iterator->raw_em_energy());
    trackster_raw_pt.push_back(trackster_iterator->raw_pt());
    trackster_raw_em_pt.push_back(trackster_iterator->raw_em_pt());
    trackster_barycenter_x.push_back(trackster_iterator->barycenter().x());
    trackster_barycenter_y.push_back(trackster_iterator->barycenter().y());
    trackster_barycenter_z.push_back(trackster_iterator->barycenter().z());
    trackster_barycenter_eta.push_back(trackster_iterator->barycenter().eta());
    trackster_barycenter_phi.push_back(trackster_iterator->barycenter().phi());
    trackster_EV1.push_back(trackster_iterator->eigenvalues()[0]);
    trackster_EV2.push_back(trackster_iterator->eigenvalues()[1]);
    trackster_EV3.push_back(trackster_iterator->eigenvalues()[2]);
    trackster_eVector0_x.push_back((trackster_iterator->eigenvectors()[0]).x());
    trackster_eVector0_y.push_back((trackster_iterator->eigenvectors()[0]).y());
    trackster_eVector0_z.push_back((trackster_iterator->eigenvectors()[0]).z());
    trackster_sigmaPCA1.push_back(trackster_iterator->sigmasPCA()[0]);
    trackster_sigmaPCA2.push_back(trackster_iterator->sigmasPCA()[1]);
    trackster_sigmaPCA3.push_back(trackster_iterator->sigmasPCA()[2]);
    std::vector<float_t> id_probs;
    for (size_t i = 0; i < 8; i++)
      id_probs.push_back(trackster_iterator->id_probabilities(i));
    trackster_id_probabilities.push_back(id_probs);
    
    // Clusters
    std::vector<uint16_t> vertices_indexes;
    std::vector<float_t> vertices_x;
    std::vector<float_t> vertices_y;
    std::vector<float_t> vertices_z;
    std::vector<float_t> vertices_energy;
    std::vector<float_t> vertices_correctedEnergy;
    std::vector<float_t> vertices_correctedEnergyUncertainty;
    for (auto idx : trackster_iterator->vertices()) {
        vertices_indexes.push_back(idx);
        auto associated_cluster = (*layer_clusters_h)[idx];
        vertices_x.push_back(associated_cluster.x());
        vertices_y.push_back(associated_cluster.y());
        vertices_z.push_back(associated_cluster.z());
        vertices_energy.push_back(associated_cluster.energy());
        vertices_correctedEnergy.push_back(associated_cluster.correctedEnergy());
        vertices_correctedEnergyUncertainty.push_back(associated_cluster.correctedEnergyUncertainty());
    }
    trackster_vertices_indexes.push_back(vertices_indexes);
    trackster_vertices_x.push_back(vertices_x);
    trackster_vertices_y.push_back(vertices_y);
    trackster_vertices_z.push_back(vertices_z);
    trackster_vertices_energy.push_back(vertices_energy);
    trackster_vertices_correctedEnergy.push_back(vertices_correctedEnergy);
    trackster_vertices_correctedEnergyUncertainty.push_back(vertices_correctedEnergyUncertainty);
    layer_cluster_seed.push_back(tracksterSeeds[t_id]);
    
    // Multiplicity
    std::vector<float_t> vertices_multiplicity;
    for (auto multiplicity : trackster_iterator->vertex_multiplicity()) {
      vertices_multiplicity.push_back(multiplicity);
    }
    trackster_vertices_multiplicity.push_back(vertices_multiplicity);
    t_id += 1;
  }

  stsSC_ntracksters_ = simTrackstersSC.size();
  nclusters_ = clusters.size();

  for (auto trackster_iterator = simTrackstersSC.begin(); trackster_iterator != simTrackstersSC.end(); ++trackster_iterator) {
      //per-trackster analysis
    stsSC_trackster_time.push_back(trackster_iterator->time());
    stsSC_trackster_timeError.push_back(trackster_iterator->timeError());
    stsSC_trackster_regressed_energy.push_back(trackster_iterator->regressed_energy());
    stsSC_trackster_raw_energy.push_back(trackster_iterator->raw_energy());
    stsSC_trackster_raw_em_energy.push_back(trackster_iterator->raw_em_energy());
    stsSC_trackster_raw_pt.push_back(trackster_iterator->raw_pt());
    stsSC_trackster_raw_em_pt.push_back(trackster_iterator->raw_em_pt());
    stsSC_trackster_barycenter_x.push_back(trackster_iterator->barycenter().x());
    stsSC_trackster_barycenter_y.push_back(trackster_iterator->barycenter().y());
    stsSC_trackster_barycenter_z.push_back(trackster_iterator->barycenter().z());
    stsSC_trackster_barycenter_eta.push_back(trackster_iterator->barycenter().eta());
    stsSC_trackster_barycenter_phi.push_back(trackster_iterator->barycenter().phi());
    stsSC_trackster_EV1.push_back(trackster_iterator->eigenvalues()[0]);
    stsSC_trackster_EV2.push_back(trackster_iterator->eigenvalues()[1]);
    stsSC_trackster_EV3.push_back(trackster_iterator->eigenvalues()[2]);
    stsSC_trackster_eVector0_x.push_back((trackster_iterator->eigenvectors()[0]).x());
    stsSC_trackster_eVector0_y.push_back((trackster_iterator->eigenvectors()[0]).y());
    stsSC_trackster_eVector0_z.push_back((trackster_iterator->eigenvectors()[0]).z());
    stsSC_trackster_sigmaPCA1.push_back(trackster_iterator->sigmasPCA()[0]);
    stsSC_trackster_sigmaPCA2.push_back(trackster_iterator->sigmasPCA()[1]);
    stsSC_trackster_sigmaPCA3.push_back(trackster_iterator->sigmasPCA()[2]);
    std::vector<float_t> id_probs;
    for (size_t i = 0; i < 8; i++)
      id_probs.push_back(trackster_iterator->id_probabilities(i));
    stsSC_trackster_id_probabilities.push_back(id_probs);
    
    // Clusters
    std::vector<uint16_t> vertices_indexes;
    std::vector<float_t> vertices_x;
    std::vector<float_t> vertices_y;
    std::vector<float_t> vertices_z;
    std::vector<float_t> vertices_energy;
    std::vector<float_t> vertices_correctedEnergy;
    std::vector<float_t> vertices_correctedEnergyUncertainty;
    for (auto idx : trackster_iterator->vertices()) {
        vertices_indexes.push_back(idx);
        auto associated_cluster = (*layer_clusters_h)[idx];
        vertices_x.push_back(associated_cluster.x());
        vertices_y.push_back(associated_cluster.y());
        vertices_z.push_back(associated_cluster.z());
        vertices_energy.push_back(associated_cluster.energy());
        vertices_correctedEnergy.push_back(associated_cluster.correctedEnergy());
        vertices_correctedEnergyUncertainty.push_back(associated_cluster.correctedEnergyUncertainty());
    }
    stsSC_trackster_vertices_indexes.push_back(vertices_indexes);
    stsSC_trackster_vertices_x.push_back(vertices_x);
    stsSC_trackster_vertices_y.push_back(vertices_y);
    stsSC_trackster_vertices_z.push_back(vertices_z);
    stsSC_trackster_vertices_energy.push_back(vertices_energy);
    stsSC_trackster_vertices_correctedEnergy.push_back(vertices_correctedEnergy);
    stsSC_trackster_vertices_correctedEnergyUncertainty.push_back(vertices_correctedEnergyUncertainty);
    
    // Multiplicity
    std::vector<float_t> vertices_multiplicity;
    for (auto multiplicity : trackster_iterator->vertex_multiplicity()) {
      vertices_multiplicity.push_back(multiplicity);
    }
    stsSC_trackster_vertices_multiplicity.push_back(vertices_multiplicity);
  }


  stsCP_ntracksters_ = simTrackstersCP.size();
  nclusters_ = clusters.size();
  for (auto trackster_iterator = simTrackstersCP.begin(); trackster_iterator != simTrackstersCP.end(); ++trackster_iterator) {
      //per-trackster analysis
    stsCP_trackster_time.push_back(trackster_iterator->time());
    stsCP_trackster_timeError.push_back(trackster_iterator->timeError());
    stsCP_trackster_regressed_energy.push_back(trackster_iterator->regressed_energy());
    stsCP_trackster_raw_energy.push_back(trackster_iterator->raw_energy());
    stsCP_trackster_raw_em_energy.push_back(trackster_iterator->raw_em_energy());
    stsCP_trackster_raw_pt.push_back(trackster_iterator->raw_pt());
    stsCP_trackster_raw_em_pt.push_back(trackster_iterator->raw_em_pt());
    stsCP_trackster_barycenter_x.push_back(trackster_iterator->barycenter().x());
    stsCP_trackster_barycenter_y.push_back(trackster_iterator->barycenter().y());
    stsCP_trackster_barycenter_z.push_back(trackster_iterator->barycenter().z());
    stsCP_trackster_barycenter_eta.push_back(trackster_iterator->barycenter().eta());
    stsCP_trackster_barycenter_phi.push_back(trackster_iterator->barycenter().phi());
    stsCP_trackster_EV1.push_back(trackster_iterator->eigenvalues()[0]);
    stsCP_trackster_EV2.push_back(trackster_iterator->eigenvalues()[1]);
    stsCP_trackster_EV3.push_back(trackster_iterator->eigenvalues()[2]);
    stsCP_trackster_eVector0_x.push_back((trackster_iterator->eigenvectors()[0]).x());
    stsCP_trackster_eVector0_y.push_back((trackster_iterator->eigenvectors()[0]).y());
    stsCP_trackster_eVector0_z.push_back((trackster_iterator->eigenvectors()[0]).z());
    stsCP_trackster_sigmaPCA1.push_back(trackster_iterator->sigmasPCA()[0]);
    stsCP_trackster_sigmaPCA2.push_back(trackster_iterator->sigmasPCA()[1]);
    stsCP_trackster_sigmaPCA3.push_back(trackster_iterator->sigmasPCA()[2]);
    std::vector<float_t> id_probs;
    for (size_t i = 0; i < 8; i++)
      id_probs.push_back(trackster_iterator->id_probabilities(i));
    stsCP_trackster_id_probabilities.push_back(id_probs);
    
    // Clusters
    std::vector<uint16_t> vertices_indexes;
    std::vector<float_t> vertices_x;
    std::vector<float_t> vertices_y;
    std::vector<float_t> vertices_z;
    std::vector<float_t> vertices_energy;
    std::vector<float_t> vertices_correctedEnergy;
    std::vector<float_t> vertices_correctedEnergyUncertainty;
    for (auto idx : trackster_iterator->vertices()) {
        vertices_indexes.push_back(idx);
        auto associated_cluster = (*layer_clusters_h)[idx];
        vertices_x.push_back(associated_cluster.x());
        vertices_y.push_back(associated_cluster.y());
        vertices_z.push_back(associated_cluster.z());
        vertices_energy.push_back(associated_cluster.energy());
        vertices_correctedEnergy.push_back(associated_cluster.correctedEnergy());
        vertices_correctedEnergyUncertainty.push_back(associated_cluster.correctedEnergyUncertainty());
    }
    stsCP_trackster_vertices_indexes.push_back(vertices_indexes);
    stsCP_trackster_vertices_x.push_back(vertices_x);
    stsCP_trackster_vertices_y.push_back(vertices_y);
    stsCP_trackster_vertices_z.push_back(vertices_z);
    stsCP_trackster_vertices_energy.push_back(vertices_energy);
    stsCP_trackster_vertices_correctedEnergy.push_back(vertices_correctedEnergy);
    stsCP_trackster_vertices_correctedEnergyUncertainty.push_back(vertices_correctedEnergyUncertainty);
    
    // Multiplicity
    std::vector<float_t> vertices_multiplicity;
    for (auto multiplicity : trackster_iterator->vertex_multiplicity()) {
      vertices_multiplicity.push_back(multiplicity);
    }
    stsCP_trackster_vertices_multiplicity.push_back(vertices_multiplicity);
  }

  
  node_linked_inners.resize(tracksters.size());
  node_linked_outers.resize(tracksters.size());
  isRootTrackster.resize(tracksters.size(),false);
  for (size_t i = 0; i < tracksters.size(); ++i) {
    const auto &node = graph.getNode((int) i);
    auto this_inners = node.getInner();
    auto this_outers = node.getOuter();
    node_linked_inners[i].insert(node_linked_inners[i].end(), this_inners.begin(), this_inners.end());
    node_linked_outers[i].insert(node_linked_outers[i].end(), this_outers.begin(), this_outers.end());
    if (node.getInner().empty()) isRootTrackster[i] = true;
  }

  int c_id = 0;
  
  for (auto cluster_iterator = clusters.begin(); cluster_iterator != clusters.end(); ++cluster_iterator) {
    auto lc_seed = cluster_iterator->seed();
    cluster_seedID.push_back(lc_seed);
    cluster_energy.push_back(cluster_iterator->energy());
    cluster_correctedEnergy.push_back(cluster_iterator->correctedEnergy());
    cluster_correctedEnergyUncertainty.push_back(cluster_iterator->correctedEnergyUncertainty());
    cluster_position_x.push_back(cluster_iterator->x());
    cluster_position_y.push_back(cluster_iterator->y());
    cluster_position_z.push_back(cluster_iterator->z());
    cluster_position_eta.push_back(cluster_iterator->eta());
    cluster_position_phi.push_back(cluster_iterator->phi());
    uint32_t number_of_hits = cluster_iterator->hitsAndFractions().size();
    cluster_number_of_hits.push_back(number_of_hits);    
    cluster_number_of_hits.push_back(number_of_hits);  
    cluster_type.push_back(ticl::returnIndex(lc_seed, rhtools_));

    cluster_timeErr.push_back(layerClustersTimes.get(c_id).second);
    cluster_time.push_back(layerClustersTimes.get(c_id).first);

    auto c_ld = 0.;
    if(layer_cluster_density[c_id] > 1e-5 && layer_cluster_density[c_id] < 3000){
      c_ld = layer_cluster_density[c_ld];
    }
    cluster_ld.push_back(c_ld);
    
    auto c_r = 0.;
    if(layer_cluster_radius[c_id] > 1e-5 && layer_cluster_radius[c_id] < 3000){
      c_r = layer_cluster_radius[c_id];
    }
    cluster_radius.push_back(c_r);
    c_id += 1;

  }
  
  tracksters_in_candidate.resize(ticlcandidates.size());
  track_in_candidate.resize(ticlcandidates.size(), -1);
  nCandidates = ticlcandidates.size();
  for (size_t i = 0; i < ticlcandidates.size(); ++i) {
    const auto& candidate = ticlcandidates[i];
    candidate_charge.push_back(candidate.charge());
    candidate_pdgId.push_back(candidate.pdgId());
    candidate_energy.push_back(candidate.energy());
    candidate_px.push_back(candidate.px());
    candidate_py.push_back(candidate.py());
    candidate_pz.push_back(candidate.pz());
    candidate_time.push_back(candidate.time());
    candidate_time_err.push_back(candidate.timeError());
    std::vector<float_t> id_probs;
    for (int i = 0; i < 8; i++) {
      ticl::Trackster::ParticleType type = static_cast<ticl::Trackster::ParticleType>(i);
      id_probs.push_back(candidate.id_probability(type));
    }
    candidate_id_probabilities.push_back(id_probs);

    auto trackster_ptrs = candidate.tracksters();
    auto track_ptr = candidate.trackPtr();
    for (auto ts_ptr : trackster_ptrs) {
      auto ts_idx = ts_ptr.get() - (edm::Ptr<ticl::Trackster>(tracksters_handle, 0)).get();
      tracksters_in_candidate[i].push_back(ts_idx);
    }
    
    if (track_ptr.isNull()) continue;
    uint16_t tk_idx = track_ptr.get() - (edm::Ptr<reco::Track>(tracks_h, 0)).get();
    track_in_candidate[i] = tk_idx;
    if(mask_tracks[tk_idx] == false){
      std::cout << "SOMETHING WRONG "<< std::endl;
    }
  }

  
  nTrackstersMerged = trackstersmerged.size();
  for (size_t i = 0; i < trackstersmerged.size(); ++i) {
    const auto& tsm = trackstersmerged[i];
    tracksters_merged_barycenter_x.push_back(tsm.barycenter().x());
    tracksters_merged_barycenter_y.push_back(tsm.barycenter().y());
    tracksters_merged_barycenter_z.push_back(tsm.barycenter().z());
    tracksters_merged_EV1.push_back(tsm.eigenvalues()[0]);
    tracksters_merged_EV2.push_back(tsm.eigenvalues()[1]);
    tracksters_merged_EV3.push_back(tsm.eigenvalues()[2]);
    tracksters_merged_eVector0_x.push_back(tsm.eigenvectors(0).x());
    tracksters_merged_eVector0_y.push_back(tsm.eigenvectors(0).y());
    tracksters_merged_eVector0_z.push_back(tsm.eigenvectors(0).z());
    tracksters_merged_sigmaPCA1.push_back(tsm.sigmasPCA()[0]);
    tracksters_merged_sigmaPCA2.push_back(tsm.sigmasPCA()[1]);
    tracksters_merged_sigmaPCA3.push_back(tsm.sigmasPCA()[2]);
    
    std::vector<float_t> id_probs;
    for (size_t i = 0; i < 8; i++)
      id_probs.push_back(tsm.id_probabilities(i));
    tracksters_merged_id_probabilities.push_back(id_probs);
  }

  
  // Tackster reco->sim associations
  trackstersCLUE3D_recoToSim_SC.resize(tracksters.size());
  trackstersCLUE3D_recoToSim_SC_score.resize(tracksters.size());
  for (size_t i = 0; i < tracksters.size(); ++i) {
    const edm::Ref<ticl::TracksterCollection> tsRef(tracksters_handle, i);
    
    // CLUE3D -> STS-SC
    const auto stsSC_iter = tsRecoSimSCMap.find(tsRef);
    if(stsSC_iter != tsRecoSimSCMap.end()) {
      const auto& stsSCassociated = stsSC_iter->val;
      for (auto &sts : stsSCassociated) {
        auto sts_id = (sts.first).get() - (edm::Ref<ticl::TracksterCollection>(simTrackstersSC_h, 0)).get();
        trackstersCLUE3D_recoToSim_SC[i].push_back(sts_id);
        trackstersCLUE3D_recoToSim_SC_score[i].push_back(sts.second);
      }
    }
  }

  
  // SimTracksters
  nsimTrackstersSC = simTrackstersSC.size();
  trackstersCLUE3D_simToReco_SC.resize(nsimTrackstersSC);
  trackstersCLUE3D_simToReco_SC_score.resize(nsimTrackstersSC);
  trackstersCLUE3D_simToReco_SC_sharedE.resize(nsimTrackstersSC);
  for (size_t i = 0; i < nsimTrackstersSC; ++i) {
    const edm::Ref<ticl::TracksterCollection> stsSCRef(simTrackstersSC_h, i);

    // STS-SC -> CLUE3D
    const auto ts_iter = tsSimToRecoSCMap.find(stsSCRef);
    if (ts_iter != tsSimToRecoSCMap.end()) {
      const auto& tsAssociated = ts_iter->val;
      for (auto &ts : tsAssociated) {
        auto ts_idx = (ts.first).get() - (edm::Ref<ticl::TracksterCollection>(tracksters_handle, i)).get();
        trackstersCLUE3D_simToReco_SC[i].push_back(ts_idx);
        trackstersCLUE3D_simToReco_SC_score[i].push_back(ts.second.second);
        trackstersCLUE3D_simToReco_SC_sharedE[i].push_back(ts.second.first);
      }
    }
  }

  
  // Tackster reco->sim associations
  trackstersCLUE3D_recoToSim_CP.resize(tracksters.size());
  trackstersCLUE3D_recoToSim_CP_score.resize(tracksters.size());
  for (size_t i = 0; i < tracksters.size(); ++i) {
    const edm::Ref<ticl::TracksterCollection> tsRef(tracksters_handle, i);
    
    // CLUE3D -> STS-CP
    const auto stsCP_iter = tsRecoSimCPMap.find(tsRef);
    if(stsCP_iter != tsRecoSimCPMap.end()) {
      const auto& stsCPassociated = stsCP_iter->val;
      for (auto &sts : stsCPassociated) {
        auto sts_id = (sts.first).get() - (edm::Ref<ticl::TracksterCollection>(simTrackstersCP_h, 0)).get();
        trackstersCLUE3D_recoToSim_CP[i].push_back(sts_id);
        trackstersCLUE3D_recoToSim_CP_score[i].push_back(sts.second);
      }
    }
  } 
  

  // SimTracksters
  nsimTrackstersCP = simTrackstersCP.size();
  trackstersCLUE3D_simToReco_CP.resize(nsimTrackstersCP);
  trackstersCLUE3D_simToReco_CP_score.resize(nsimTrackstersCP);
  trackstersCLUE3D_simToReco_CP_sharedE.resize(nsimTrackstersCP);
  for (size_t i = 0; i < nsimTrackstersCP; ++i) {
    const edm::Ref<ticl::TracksterCollection> stsCPRef(simTrackstersCP_h, i);

    // STS-CP -> CLUE3D
    const auto ts_iter = tsSimToRecoCPMap.find(stsCPRef);
    if (ts_iter != tsSimToRecoCPMap.end()) {
      const auto& tsAssociated = ts_iter->val;
      for (auto &ts : tsAssociated) {
        auto ts_idx = (ts.first).get() - (edm::Ref<ticl::TracksterCollection>(tracksters_handle, i)).get();
        trackstersCLUE3D_simToReco_CP[i].push_back(ts_idx);
        trackstersCLUE3D_simToReco_CP_score[i].push_back(ts.second.second);
        trackstersCLUE3D_simToReco_CP_sharedE[i].push_back(ts.second.first);
      }
    }
  }
  
  // Tackster reco->sim associations
  MergeTracksters_recoToSim_SC.resize(trackstersmerged.size());
  MergeTracksters_recoToSim_SC_score.resize(trackstersmerged.size());
  for (size_t i = 0; i < trackstersmerged.size(); ++i) {
    const edm::Ref<ticl::TracksterCollection> tsRef(tracksters_merged_h, i);
    
    // CLUE3D -> STS-SC
    const auto stsSC_iter = MergetsRecoSimSCMap.find(tsRef);
    if(stsSC_iter != MergetsRecoSimSCMap.end()) {
      const auto& stsSCassociated = stsSC_iter->val;
      for (auto &sts : stsSCassociated) {
        auto sts_id = (sts.first).get() - (edm::Ref<ticl::TracksterCollection>(simTrackstersSC_h, 0)).get();
        MergeTracksters_recoToSim_SC[i].push_back(sts_id);
        MergeTracksters_recoToSim_SC_score[i].push_back(sts.second);
      }
    }
  }
  
  // SimTracksters
  nsimTrackstersSC = simTrackstersSC.size();
  MergeTracksters_simToReco_SC.resize(nsimTrackstersSC);
  MergeTracksters_simToReco_SC_score.resize(nsimTrackstersSC);
  MergeTracksters_simToReco_SC_sharedE.resize(nsimTrackstersSC);
  for (size_t i = 0; i < nsimTrackstersSC; ++i) {
    const edm::Ref<ticl::TracksterCollection> stsSCRef(simTrackstersSC_h, i);

    // STS-SC -> CLUE3D
    const auto ts_iter = MergetsSimToRecoSCMap.find(stsSCRef);
    if (ts_iter != MergetsSimToRecoSCMap.end()) {
      const auto& tsAssociated = ts_iter->val;
      for (auto &ts : tsAssociated) {
        auto ts_idx = (ts.first).get() - (edm::Ref<ticl::TracksterCollection>(tracksters_merged_h, i)).get();
        MergeTracksters_simToReco_SC[i].push_back(ts_idx);
        MergeTracksters_simToReco_SC_score[i].push_back(ts.second.second);
        MergeTracksters_simToReco_SC_sharedE[i].push_back(ts.second.first);
      }
    }
  }

  
  // Tackster reco->sim associations
  MergeTracksters_recoToSim_CP.resize(trackstersmerged.size());
  MergeTracksters_recoToSim_CP_score.resize(trackstersmerged.size());
  for (size_t i = 0; i < trackstersmerged.size(); ++i) {
    const edm::Ref<ticl::TracksterCollection> tsRef(tracksters_merged_h, i);
    
    // CLUE3D -> STS-CP
    const auto stsCP_iter = MergetsRecoSimCPMap.find(tsRef);
    if(stsCP_iter != MergetsRecoSimCPMap.end()) {
      const auto& stsCPassociated = stsCP_iter->val;
      for (auto &sts : stsCPassociated) {
        auto sts_id = (sts.first).get() - (edm::Ref<ticl::TracksterCollection>(simTrackstersCP_h, 0)).get();
        MergeTracksters_recoToSim_CP[i].push_back(sts_id);
        MergeTracksters_recoToSim_CP_score[i].push_back(sts.second);
      }
    }
  }


  
  // SimTracksters
  nsimTrackstersCP = simTrackstersCP.size();
  MergeTracksters_simToReco_CP.resize(nsimTrackstersCP);
  MergeTracksters_simToReco_CP_score.resize(nsimTrackstersCP);
  MergeTracksters_simToReco_CP_sharedE.resize(nsimTrackstersCP);
  for (size_t i = 0; i < nsimTrackstersCP; ++i) {
    const edm::Ref<ticl::TracksterCollection> stsCPRef(simTrackstersCP_h, i);
  
    // STS-CP -> CLUE3D
    const auto ts_iter = MergetsSimToRecoCPMap.find(stsCPRef);
    if (ts_iter != MergetsSimToRecoCPMap.end()) {
      const auto& tsAssociated = ts_iter->val;
      for (auto &ts : tsAssociated) {
        auto ts_idx = (ts.first).get() - (edm::Ref<ticl::TracksterCollection>(tracksters_merged_h, i)).get();
        MergeTracksters_simToReco_CP[i].push_back(ts_idx);
        MergeTracksters_simToReco_CP_score[i].push_back(ts.second.second);
        MergeTracksters_simToReco_CP_sharedE[i].push_back(ts.second.first);
      }
    }
  }

  
  //Tracks 
  for(size_t i = 0; i < tracks.size(); i++){
    if(mask_tracks[i] == true){
      reco::TrackRef trackref = reco::TrackRef(tracks_h, i);
      
      track_ev.push_back(event_index);
      track_id.push_back(i);  

      track_hgcal_x.push_back(tracks_propagated_x[i]);

      track_hgcal_y.push_back(tracks_propagated_y[i]);

      track_hgcal_z.push_back(tracks_propagated_z[i]);

      track_hgcal_eta.push_back(tracks_propagated_eta[i]);

      track_hgcal_phi.push_back(tracks_propagated_phi[i]);
      
      track_hgcal_px.push_back(tracks_propagated_px[i]);
      
      track_hgcal_py.push_back(tracks_propagated_py[i]);
      
      track_hgcal_pz.push_back(tracks_propagated_pz[i]);

      
      track_pt.push_back(tracks[i].pt());
      track_charge.push_back(tracks[i].charge());
      
      track_time.push_back(trackTime[trackref]);      
      track_time_quality.push_back(trackTimeQual[trackref]);
      track_time_err.push_back(trackTimeErr[trackref]);
      track_nhits.push_back(tracks[i].recHitsSize());
    }
  }

  
  trackster_tree_->Fill();
  cluster_tree_->Fill();
  graph_tree_->Fill();
  candidate_tree_->Fill();
  tracksters_merged_tree_->Fill();
  associations_tree_->Fill();
  simtrackstersSC_tree_->Fill();
  simtrackstersCP_tree_->Fill();
  tracks_tree_->Fill();
}

void Ntupler::endJob() {}

void Ntupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackstersclue3d", edm::InputTag("ticlTrackstersCLUE3DHigh"));
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<edm::InputTag>("layer_clustersTime", edm::InputTag("hgcalLayerClusters", "timeLayerCluster"));
  desc.add<edm::InputTag>("ticlgraph", edm::InputTag("ticlGraph"));
  desc.add<edm::InputTag>("ticlcandidates", edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("masked_tracks", edm::InputTag("ticlTrackstersMerge", "maskTracks"));
  desc.add<edm::InputTag>("hgcaltracks_x", edm::InputTag("ticlTrackstersMerge", "hgcaltracksX"));
  desc.add<edm::InputTag>("hgcaltracks_y", edm::InputTag("ticlTrackstersMerge", "hgcaltracksY"));
  desc.add<edm::InputTag>("hgcaltracks_z", edm::InputTag("ticlTrackstersMerge", "hgcaltracksZ"));
  desc.add<edm::InputTag>("hgcaltracks_eta", edm::InputTag("ticlTrackstersMerge", "hgcaltracksEta"));
  desc.add<edm::InputTag>("hgcaltracks_phi", edm::InputTag("ticlTrackstersMerge", "hgcaltracksPhi"));
  desc.add<edm::InputTag>("hgcaltracks_px", edm::InputTag("ticlTrackstersMerge", "hgcaltracksPx"));
  desc.add<edm::InputTag>("hgcaltracks_py", edm::InputTag("ticlTrackstersMerge", "hgcaltracksPx"));
  desc.add<edm::InputTag>("hgcaltracks_pz", edm::InputTag("ticlTrackstersMerge", "hgcaltracksPx"));
  desc.add<edm::InputTag>("tracksTime", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("tracksTimeQual", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<edm::InputTag>("tracksTimeErr", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("trackstersmerged", edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("layerClustersLocalDensity", edm::InputTag("ticlTrackstersCLUE3DHigh", "layerClustersLocalDensity"));
  desc.add<edm::InputTag>("layerClustersRadius", edm::InputTag("ticlTrackstersCLUE3DHigh", "layerClustersRadius"));
  desc.add<edm::InputTag>("tracksterSeeds", edm::InputTag("ticlTrackstersCLUE3DHigh", "tracksterSeeds"));
  desc.add<edm::InputTag>("simtrackstersSC", edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>("simtrackstersCP", edm::InputTag("ticlSimTracksters", "fromCPs"));
  desc.add<edm::InputTag>("recoToSimAssociatorSC", edm::InputTag("tracksterSimTracksterAssociationPRbyCLUE3D"));
  desc.add<edm::InputTag>("simToRecoAssociatorSC", edm::InputTag("tracksterSimTracksterAssociationPRbyCLUE3D"));
  desc.add<edm::InputTag>("recoToSimAssociatorCP", edm::InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3D"));
  desc.add<edm::InputTag>("simToRecoAssociatorCP", edm::InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3D"));
  desc.add<edm::InputTag>("MergerecoToSimAssociatorSC", edm::InputTag("tracksterSimTracksterAssociationPR"));
  desc.add<edm::InputTag>("MergesimToRecoAssociatorSC", edm::InputTag("tracksterSimTracksterAssociationPR"));
  desc.add<edm::InputTag>("MergerecoToSimAssociatorCP", edm::InputTag("tracksterSimTracksterAssociationLinking"));
  desc.add<edm::InputTag>("MergesimToRecoAssociatorCP", edm::InputTag("tracksterSimTracksterAssociationLinking"));
  descriptions.add("ticlNtuplizer", desc);
}

DEFINE_FWK_MODULE(Ntupler);
