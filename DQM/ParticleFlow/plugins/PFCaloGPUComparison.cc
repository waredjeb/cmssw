#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitHostCollection.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterHostCollection.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#ifdef PFLOW_DEBUG
#define LOGVERB(x) edm::LogVerbatim(x)
#else
#define LOGVERB(x) LogTrace(x)
#endif

class PFCaloGPUComparison : public DQMEDAnalyzer {
public:
  PFCaloGPUComparison(edm::ParameterSet const& conf);
  ~PFCaloGPUComparison() override = default;
  void analyze(edm::Event const& e, edm::EventSetup const& c) override;
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

private:
  const edm::EDGetTokenT<reco::PFClusterHostCollection> pfClusterTok_ref_;
  const edm::EDGetTokenT<reco::PFClusterHostCollection> pfClusterTok_target_;
  const edm::EDGetTokenT<reco::PFRecHitHostCollection> pfRecHitsTok_ref_;
  const edm::EDGetTokenT<reco::PFRecHitHostCollection> pfRecHitsTok_target_;
//  edm::EDGetTokenT<reco::PFClusterCollection> pfClusterTok_ref_;
//  edm::EDGetTokenT<reco::PFClusterCollection> pfClusterTok_target_;

  MonitorElement* pfCluster_Multiplicity_GPUvsCPU_;
  MonitorElement* pfCluster_Energy_GPUvsCPU_;
  MonitorElement* pfCluster_RecHitMultiplicity_GPUvsCPU_;
  MonitorElement* pfCluster_Depth_GPUvsCPU_;
  MonitorElement* pfCluster_x_GPUvsCPU_;
  MonitorElement* pfCluster_y_GPUvsCPU_;
  MonitorElement* pfCluster_z_GPUvsCPU_;
  MonitorElement* pfCluster_DuplicateMatches_GPUvsCPU_;

  std::string pfCaloGPUCompDir_;
  std::string refStr_;
  std::string targetStr_;
};

PFCaloGPUComparison::PFCaloGPUComparison(const edm::ParameterSet& conf)
    : pfClusterTok_ref_{consumes<reco::PFClusterHostCollection>(
          conf.getUntrackedParameter<edm::InputTag>("pfClusterToken_ref"))},
      pfClusterTok_target_{consumes<reco::PFClusterHostCollection>(
          conf.getUntrackedParameter<edm::InputTag>("pfClusterToken_target"))},
     pfRecHitsTok_ref_{consumes<reco::PFRecHitHostCollection>(
          conf.getUntrackedParameter<edm::InputTag>("pfRecHitsToken_ref"))},
     pfRecHitsTok_target_{consumes<reco::PFRecHitHostCollection>(
          conf.getUntrackedParameter<edm::InputTag>("pfRecHitsToken_target"))},
      pfCaloGPUCompDir_{conf.getUntrackedParameter<std::string>("pfCaloGPUCompDir")},
      refStr_{conf.getUntrackedParameter<std::string>("GPU")},
      targetStr_{conf.getUntrackedParameter<std::string>("CPU")}
{}

void PFCaloGPUComparison::bookHistograms(DQMStore::IBooker& ibooker,
                                             edm::Run const& irun,
                                             edm::EventSetup const& isetup) {
  std::string histo;

  ibooker.setCurrentFolder("ParticleFlow/" + pfCaloGPUCompDir_);
  
  std::string refVsTargetStr = refStr_ + targetStr_;
  histo = "pfCluster_Multiplicity_"+refVsTargetStr;
  pfCluster_Multiplicity_GPUvsCPU_ = ibooker.book2D(histo, histo, 100, 0, 2000, 100, 0, 2000);

  histo = "pfCluster_Energy_"+refVsTargetStr;
  pfCluster_Energy_GPUvsCPU_ = ibooker.book2D(histo, histo, 100, 0, 500, 100, 0, 500);

  histo = "pfCluster_RecHitMultiplicity_"+refVsTargetStr;
  pfCluster_RecHitMultiplicity_GPUvsCPU_ = ibooker.book2D(histo, histo, 100, 0, 100, 100, 0, 100);

  histo = "pfCluster_Depth_"+refVsTargetStr;
  pfCluster_Depth_GPUvsCPU_ = ibooker.book2D(histo, histo, 10, 0, 10, 10, 0, 10);

  histo = "pfCluster_x_"+refVsTargetStr;
  pfCluster_x_GPUvsCPU_ = ibooker.book2D(histo, histo, 100, -200, 200, 100, -200, 200);

  histo = "pfCluster_y_"+refVsTargetStr;
  pfCluster_y_GPUvsCPU_ = ibooker.book2D(histo, histo, 100, -200, 200, 100, -200, 200);

  histo = "pfCluster_z_"+refVsTargetStr;
  pfCluster_z_GPUvsCPU_ = ibooker.book2D(histo, histo, 100, -200, 200, 100, -200, 200);

  histo = "pfCluster_DuplicateMatches_"+refVsTargetStr;
  pfCluster_DuplicateMatches_GPUvsCPU_ = ibooker.book1D(histo, histo, 100, 0., 1000);
}
void PFCaloGPUComparison::analyze(edm::Event const& event, edm::EventSetup const& c) {

  auto const& pfClusterSoA_ref = event.get(pfClusterTok_ref_).const_view(); 
  auto const& pfClusterSoA_target  = event.get(pfClusterTok_target_).const_view(); 
  auto const& rechitsHandle_ref = event.get(pfRecHitsTok_ref_).const_view();
  auto const& rechitsHandle_target = event.get(pfRecHitsTok_target_).const_view();
  // Compare per-event PF cluster multiplicity
  //
  auto size_ref = pfClusterSoA_ref.nSeeds();
  auto size_target = pfClusterSoA_target.nSeeds();
  if(size_ref != size_target)
    LOGVERB("PFCaloGPUComparison") << " PFCluster multiplicity " << size_ref << " "
                                       << size_target;
  pfCluster_Multiplicity_GPUvsCPU_->Fill((float)size_ref, (float)size_target);
  std::vector<int> matched_idx;
  matched_idx.reserve(size_ref);
  for(int i = 0; i < size_ref; ++i){
    bool matched = false;
    auto rhSeed_ref = pfClusterSoA_ref[i].seedRHIdx();
    auto rh_ref = rechitsHandle_ref[rhSeed_ref].detId();
    for(int j = 0; j < size_target; ++j){
      auto rhSeed_target= pfClusterSoA_target[j].seedRHIdx();
      auto rh_target = rechitsHandle_target[rhSeed_target].detId();
      if(rh_ref == rh_target){
        if(!matched){
          matched = true;
          matched_idx.push_back(j);
        }
        else{
          edm::LogWarning("PFCaloGPUComparison") << "Found duplicate match";
          pfCluster_DuplicateMatches_GPUvsCPU_->Fill((int)j);
        }
      }
    }
    if(!matched)
      matched_idx.push_back(-1);
  }
  // Plot matching PF cluster variables
  for (int i = 0; i < size_ref; ++i) {
    if (matched_idx[i] >= 0) {
      auto j = matched_idx[i];
      if(j > size_target){
        std::cout << " WTF " << std::endl;
      }
      if(i > size_ref ){
        std::cout << " WTF2 " << std::endl;
      }
      int ref_energy_bin = pfCluster_Energy_GPUvsCPU_->getTH2F()->GetXaxis()->FindBin(pfClusterSoA_ref[i].energy());
      int target_energy_bin =
          pfCluster_Energy_GPUvsCPU_->getTH2F()->GetXaxis()->FindBin(pfClusterSoA_target[j].energy());
      if (ref_energy_bin != target_energy_bin){
        edm::LogPrint("PFCaloGPUComparison")
            << "Off-diagonal energy bin entries: " << pfClusterSoA_ref[i].energy() << " "
            << pfClusterSoA_ref[i].x() << " " << pfClusterSoA_ref[i].y() << " " << pfClusterSoA_ref[i].z() << " "
            << pfClusterSoA_target[j].energy() << " " 
            << pfClusterSoA_target[j].x() << " " << pfClusterSoA_target[j].y() << " " << pfClusterSoA_target[j].z(); 
      }
      else{
      std::cout << pfClusterSoA_ref[i].energy() << " " << pfClusterSoA_target[j].energy() << " " << ref_energy_bin << " " << target_energy_bin <<  std::endl;
      pfCluster_Energy_GPUvsCPU_->Fill(pfClusterSoA_ref[i].energy(), pfClusterSoA_target[j].energy());
      pfCluster_x_GPUvsCPU_->Fill(pfClusterSoA_ref[i].x(), pfClusterSoA_target[j].x());
      pfCluster_y_GPUvsCPU_->Fill(pfClusterSoA_ref[i].y(), pfClusterSoA_target[j].y());
      pfCluster_z_GPUvsCPU_->Fill(pfClusterSoA_ref[i].z(), pfClusterSoA_target[j].z());
      pfCluster_Depth_GPUvsCPU_->Fill(pfClusterSoA_ref[i].depth(), pfClusterSoA_target[j].depth());
      pfCluster_RecHitMultiplicity_GPUvsCPU_->Fill((float)pfClusterSoA_ref[i].rhfracSize(),
                                                   (float)pfClusterSoA_target[j].rhfracSize());
      }
    }
  }

//  if (pfClusters_ref->size() != pfClusters_target->size())
//    LOGVERB("PFCaloGPUComparison") << " PFCluster multiplicity " << pfClusters_ref->size() << " "
//                                       << pfClusters_target->size();
//  pfCluster_Multiplicity_GPUvsCPU_->Fill((float)pfClusters_ref->size(), (float)pfClusters_target->size());
//
//  //
//  // Find matching PF cluster pairs
//  std::vector<int> matched_idx;
//  matched_idx.reserve(pfClusters_ref->size());
//  for (unsigned i = 0; i < pfClusters_ref->size(); ++i) {
//    bool matched = false;
//    for (unsigned j = 0; j < pfClusters_target->size(); ++j) {
//      if (pfClusters_ref->at(i).seed() == pfClusters_target->at(j).seed()) {
//        if (!matched) {
//          matched = true;
//          matched_idx.push_back((int)j);
//        } else {
//          edm::LogWarning("PFCaloGPUComparison") << "Found duplicate match";
//          pfCluster_DuplicateMatches_GPUvsCPU_->Fill((int)j);
//        }
//      }
//    }
//    if (!matched)
//      matched_idx.push_back(-1);  // if you don't find a match, put a dummy number
//  }
//
//  //
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFCaloGPUComparison);
