// -*- C++ -*-
//

#include "DQM/HcalCommon/interface/DQTask.h"
#include "DQM/HcalCommon/interface/Utilities.h"
#include "DQM/HcalCommon/interface/HashFilter.h"
#include "DQM/HcalCommon/interface/Container1D.h"
#include "DQM/HcalCommon/interface/Container2D.h"
#include "DQM/HcalCommon/interface/ContainerProf1D.h"
#include "DQM/HcalCommon/interface/ContainerProf2D.h"
#include "DQM/HcalCommon/interface/ContainerSingle1D.h"
#include "DQM/HcalCommon/interface/ContainerSingle2D.h"
#include "DQM/HcalCommon/interface/ContainerSingleProf2D.h"
#include "DQM/HcalCommon/interface/ElectronicsMap.h"
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

#include <cmath>
#ifdef PFLOW_DEBUG
#define LOGVERB(x) edm::LogVerbatim(x)
#else
#define LOGVERB(x) LogTrace(x)
#endif

using namespace hcaldqm;
using namespace hcaldqm::constants;
using namespace hcaldqm::filter;

class PFHcalGPUComparisonTask : public hcaldqm::DQTask {
public:
  PFHcalGPUComparisonTask(edm::ParameterSet const&);
  ~PFHcalGPUComparisonTask() override = default;

  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  std::shared_ptr<hcaldqm::Cache> globalBeginLuminosityBlock(edm::LuminosityBlock const&,
                                                             edm::EventSetup const&) const override;
  void globalEndLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void _process(edm::Event const&, edm::EventSetup const&) override;
  void _resetMonitors(hcaldqm::UpdateFreq) override;

  edm::EDGetTokenT<reco::PFClusterCollection> pfClusterTok_ref_;
  edm::EDGetTokenT<reco::PFClusterCollection> pfClusterTok_target_;

  MonitorElement* pfCluster_Multiplicity_GPUvsCPU_;
  MonitorElement* pfCluster_Energy_GPUvsCPU_;
  MonitorElement* pfCluster_RecHitMultiplicity_GPUvsCPU_;
  MonitorElement* pfCluster_Layer_GPUvsCPU_;
  MonitorElement* pfCluster_Depth_GPUvsCPU_;
  MonitorElement* pfCluster_Eta_GPUvsCPU_;
  MonitorElement* pfCluster_Phi_GPUvsCPU_;
  MonitorElement* pfCluster_DuplicateMatches_GPUvsCPU_;

  std::string pfCaloGPUCompDir_;
};

PFHcalGPUComparisonTask::PFHcalGPUComparisonTask(edm::ParameterSet const& conf)
    : DQTask(conf),
      pfClusterTok_ref_{consumes<reco::PFClusterCollection>(
          conf.getUntrackedParameter<edm::InputTag>("pfClusterToken_ref"))},
      pfClusterTok_target_{
          consumes<reco::PFClusterCollection>(conf.getUntrackedParameter<edm::InputTag>("pfClusterToken_target"))},
      pfCaloGPUCompDir_{conf.getUntrackedParameter<std::string>("name")} {}

/* virtual */ void PFHcalGPUComparisonTask::bookHistograms(DQMStore::IBooker& ibooker,
                                                         edm::Run const& r,
                                                         edm::EventSetup const& es) {
  DQTask::bookHistograms(ibooker, r, es);

  //	GET WHAT YOU NEED
//  edm::ESHandle<HcalDbService> dbs = es.getHandle(hcalDbServiceToken_);
//  _emap = dbs->getHcalMapping();

  //	Book monitoring elements
  const char* histo;

  ibooker.setCurrentFolder("ParticleFlow/" + pfCaloGPUCompDir_);

  histo = "pfCluster_Multiplicity_GPUvsCPU";
  pfCluster_Multiplicity_GPUvsCPU_ = ibooker.book2D(histo, histo, 100, 0, 2000, 100, 0, 2000);

  histo = "pfCluster_Energy_GPUvsCPU";
  pfCluster_Energy_GPUvsCPU_ = ibooker.book2D(histo, histo, 100, 0, 500, 100, 0, 500);

  histo = "pfCluster_RecHitMultiplicity_GPUvsCPU";
  pfCluster_RecHitMultiplicity_GPUvsCPU_ = ibooker.book2D(histo, histo, 100, 0, 100, 100, 0, 100);

  histo = "pfCluster_Layer_GPUvsCPU";
  pfCluster_Layer_GPUvsCPU_ = ibooker.book2D(histo, histo, 100, 0, 100, 100, 0, 100);

  histo = "pfCluster_Depth_GPUvsCPU";
  pfCluster_Depth_GPUvsCPU_ = ibooker.book2D(histo, histo, 100, 0, 100, 100, 0, 100);

  histo = "pfCluster_Eta_GPUvsCPU";
  pfCluster_Eta_GPUvsCPU_ = ibooker.book2D(histo, histo, 100, 0, 100, 100, 0, 100);

  histo = "pfCluster_Phi_GPUvsCPU";
  pfCluster_Phi_GPUvsCPU_ = ibooker.book2D(histo, histo, 100, 0, 100, 100, 0, 100);

  histo = "pfCluster_DuplicateMatches_GPUvsCPU";
  pfCluster_DuplicateMatches_GPUvsCPU_ = ibooker.book1D(histo, histo, 100, 0., 1000);
}

/* virtual */ void PFHcalGPUComparisonTask::_resetMonitors(hcaldqm::UpdateFreq uf) { DQTask::_resetMonitors(uf); }

/* virtual */ void PFHcalGPUComparisonTask::_process(edm::Event const& event, edm::EventSetup const&) {
  edm::Handle<reco::PFClusterCollection> pfClusters_ref;
  event.getByToken(pfClusterTok_ref_, pfClusters_ref);

  edm::Handle<reco::PFClusterCollection> pfClusters_target;
  event.getByToken(pfClusterTok_target_, pfClusters_target);

  auto lumiCache = luminosityBlockCache(event.getLuminosityBlock().index());
  _currentLS = lumiCache->currentLS;
  // Compare per-event PF cluster multiplicity

  if (pfClusters_ref->size() != pfClusters_target->size())
    LOGVERB("PFCaloGPUComparisonTask") << " PFCluster multiplicity " << pfClusters_ref->size() << " "
                                       << pfClusters_target->size();
  pfCluster_Multiplicity_GPUvsCPU_->Fill((float)pfClusters_ref->size(), (float)pfClusters_target->size());

  //
  // Find matching PF cluster pairs
  std::vector<int> matched_idx;
  matched_idx.reserve(pfClusters_ref->size());
  for (unsigned i = 0; i < pfClusters_ref->size(); ++i) {
    bool matched = false;
    for (unsigned j = 0; j < pfClusters_target->size(); ++j) {
      if (pfClusters_ref->at(i).seed() == pfClusters_target->at(j).seed()) {
        if (!matched) {
          matched = true;
          matched_idx.push_back((int)j);
        } else {
          edm::LogWarning("PFCaloGPUComparisonTask") << "Found duplicate match";
          pfCluster_DuplicateMatches_GPUvsCPU_->Fill((int)j);
        }
      }
    }
    if (!matched)
      matched_idx.push_back(-1);  // if you don't find a match, put a dummy number
  }

  //
  // Plot matching PF cluster variables
  for (unsigned i = 0; i < pfClusters_ref->size(); ++i) {
    if (matched_idx[i] >= 0) {
      unsigned int j = matched_idx[i];
      int ref_energy_bin = pfCluster_Energy_GPUvsCPU_->getTH2F()->GetXaxis()->FindBin(pfClusters_ref->at(i).energy());
      int target_energy_bin =
          pfCluster_Energy_GPUvsCPU_->getTH2F()->GetXaxis()->FindBin(pfClusters_target->at(j).energy());
      if (ref_energy_bin != target_energy_bin)
        edm::LogPrint("PFCaloGPUComparisonTask")
            << "Off-diagonal energy bin entries: " << pfClusters_ref->at(i).energy() << " "
            << pfClusters_ref->at(i).eta() << " " << pfClusters_ref->at(i).phi() << " "
            << pfClusters_target->at(j).energy() << " " << pfClusters_target->at(j).eta() << " "
            << pfClusters_target->at(j).phi() << std::endl;
      pfCluster_Energy_GPUvsCPU_->Fill(pfClusters_ref->at(i).energy(), pfClusters_target->at(j).energy());
      pfCluster_Layer_GPUvsCPU_->Fill(pfClusters_ref->at(i).layer(), pfClusters_target->at(j).layer());
      pfCluster_Eta_GPUvsCPU_->Fill(pfClusters_ref->at(i).eta(), pfClusters_target->at(j).eta());
      pfCluster_Phi_GPUvsCPU_->Fill(pfClusters_ref->at(i).phi(), pfClusters_target->at(j).phi());
      pfCluster_Depth_GPUvsCPU_->Fill(pfClusters_ref->at(i).depth(), pfClusters_target->at(j).depth());
      pfCluster_RecHitMultiplicity_GPUvsCPU_->Fill((float)pfClusters_ref->at(i).recHitFractions().size(),
                                                   (float)pfClusters_target->at(j).recHitFractions().size());
    }
  }

  //std::map<HcalDetId, double> mRecHitEnergy;

  //for (HBHERecHitCollection::const_iterator it = chbhe_ref->begin(); it != chbhe_ref->end(); ++it) {
  //  double energy = it->energy();

  //  //	Explicit check on the DetIds present in the Collection
  //  HcalDetId did = it->id();

  //  if (mRecHitEnergy.find(did) == mRecHitEnergy.end())
  //    mRecHitEnergy.insert(std::make_pair(did, energy));
  //  else
  //    edm::LogError("PFHcalGPUComparisonTask") << "Duplicate Rechit from the same HcalDetId";
  //  ;
  //}

  //for (HBHERecHitCollection::const_iterator it = chbhe_target->begin(); it != chbhe_target->end(); ++it) {
  //  double energy = it->energy();
  //  HcalDetId did = it->id();

  //  if (mRecHitEnergy.find(did) != mRecHitEnergy.end()) {
  //    energyGPUvsCPU_subdet_.fill(did, mRecHitEnergy[did], energy);

  //    if (mRecHitEnergy[did] != 0.) {
  //      energyDiffGPUCPU_subdet_.fill(did, (energy - mRecHitEnergy[did]) / mRecHitEnergy[did]);
  //      if (energy > 0.1)
  //        energyDiffGPUCPU_depth_.fill(did, (energy - mRecHitEnergy[did]) / mRecHitEnergy[did]);
  //    } else if (mRecHitEnergy[did] == 0. && energy == 0.) {
  //      energyDiffGPUCPU_subdet_.fill(did, 0.);
  //      if (energy > 0.1)
  //        energyDiffGPUCPU_depth_.fill(did, 0.);
  //    } else {
  //      energyDiffGPUCPU_subdet_.fill(did, -1.);
  //      if (energy > 0.1)
  //        energyDiffGPUCPU_depth_.fill(did, -1.);
  //    }

  //    mRecHitEnergy.erase(did);
  //  } else {
  //    if (energy > 2.)
  //      edm::LogError("PFHcalGPUComparisonTask")
  //          << "Energetic GPU Rechit exist, but not reconstructed by CPU. DetId = " << did;
  //  }
  //}
  //if (!mRecHitEnergy.empty()) {
  //  for (auto const& rhpair : mRecHitEnergy) {
  //    if (rhpair.second > 2.)
  //      edm::LogError("PFHcalGPUComparisonTask")
  //          << "Energetic CPU Rechit exist, but not reconstructed by GPU. DetId = " << rhpair.first;
  //  }
  //}
}

std::shared_ptr<hcaldqm::Cache> PFHcalGPUComparisonTask::globalBeginLuminosityBlock(edm::LuminosityBlock const& lb,
                                                                                  edm::EventSetup const& es) const {
  return DQTask::globalBeginLuminosityBlock(lb, es);
}

/* virtual */ void PFHcalGPUComparisonTask::globalEndLuminosityBlock(edm::LuminosityBlock const& lb,
                                                                   edm::EventSetup const& es) {
  if (_ptype != fOnline)
    return;

  auto lumiCache = luminosityBlockCache(lb.index());
  _currentLS = lumiCache->currentLS;

  //	in the end always do the DQTask::endLumi
  DQTask::globalEndLuminosityBlock(lb, es);
}

void PFHcalGPUComparisonTask::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<std::string>("name", "pfCaloGPUCompDir");
  desc.addUntracked<edm::InputTag>("pfClusterToken_ref", edm::InputTag("hltParticleFlowClusterHCALSerialSync"));
  desc.addUntracked<edm::InputTag>("pfClusterToken_target", edm::InputTag("hltParticleFlowClusterHCAL"));
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(PFHcalGPUComparisonTask);
