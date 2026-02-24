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
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
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

class HGCalLayerClustersGPUComparisonTask : public hcaldqm::DQTask {
public:
  HGCalLayerClustersGPUComparisonTask(edm::ParameterSet const&);
  ~HGCalLayerClustersGPUComparisonTask() override = default;

  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  std::shared_ptr<hcaldqm::Cache> globalBeginLuminosityBlock(edm::LuminosityBlock const&,
                                                             edm::EventSetup const&) const override;
  void globalEndLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void _process(edm::Event const&, edm::EventSetup const&) override;
  void _resetMonitors(hcaldqm::UpdateFreq) override;

  edm::EDGetTokenT<std::vector<reco::CaloCluster>> layerClusterTok_ref_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> layerClusterTok_target_;

  MonitorElement* layerCluster_Multiplicity_HostvsDevice_;
  MonitorElement* layerCluster_Energy_HostvsDevice_;
  MonitorElement* layerCluster_RecHitMultiplicity_HostvsDevice_;
  MonitorElement* layerCluster_Layer_HostvsDevice_;
  MonitorElement* layerCluster_Depth_HostvsDevice_;
  MonitorElement* layerCluster_Eta_HostvsDevice_;
  MonitorElement* layerCluster_Phi_HostvsDevice_;
  MonitorElement* layerCluster_DuplicateMatches_HostvsDevice_;

  MonitorElement* layerCluster_Multiplicity_Diff_HostvsDevice_;
  MonitorElement* layerCluster_Energy_Diff_HostvsDevice_;
  MonitorElement* layerCluster_RecHitMultiplicity_Diff_HostvsDevice_;
  MonitorElement* layerCluster_Layer_Diff_HostvsDevice_;
  MonitorElement* layerCluster_Depth_Diff_HostvsDevice_;
  MonitorElement* layerCluster_Eta_Diff_HostvsDevice_;
  MonitorElement* layerCluster_Phi_Diff_HostvsDevice_;

  std::string subsystemDir_;
  std::string pfCaloGPUCompDir_;
};

HGCalLayerClustersGPUComparisonTask::HGCalLayerClustersGPUComparisonTask(edm::ParameterSet const& conf)
    : DQTask(conf),
      layerClusterTok_ref_{
          consumes<std::vector<reco::CaloCluster>>(conf.getUntrackedParameter<edm::InputTag>("layerClusterToken_ref"))},
      layerClusterTok_target_{
          consumes<std::vector<reco::CaloCluster>>(conf.getUntrackedParameter<edm::InputTag>("layerClusterToken_target"))},
      subsystemDir_{conf.getUntrackedParameter<std::string>("subsystem")},
      pfCaloGPUCompDir_{conf.getUntrackedParameter<std::string>("name")} {}

void HGCalLayerClustersGPUComparisonTask::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const& r, edm::EventSetup const& es) {
  _subsystem = subsystemDir_;
  ibooker.setCurrentFolder(pfCaloGPUCompDir_);
  DQTask::bookHistograms(ibooker, r, es);
  //	Book monitoring elements
  const char* histo;

  histo = "layerCluster_Multiplicity_HostvsDevice";
  const char* histoAxis = "layerCluster_Multiplicity_HostvsDevice;Multiplicity Device;Multiplicity Device";
  layerCluster_Multiplicity_HostvsDevice_ = ibooker.book2I(histo, histoAxis, 1000, 0, 1000, 1000, 0, 1000);

  histo = "layerCluster_Energy_HostvsDevice";
  histoAxis = "layerCluster_Energy_HostvsDevice;Energy Host [GeV];Energy Device [GeV]";
  layerCluster_Energy_HostvsDevice_ = ibooker.book2D(histo, histoAxis, 500, 0, 500, 500, 0, 500);

  histo = "layerCluster_RecHitMultiplicity_HostvsDevice";
  histoAxis = "layerCluster_RecHitMultiplicity_HostvsDevice;RecHit Multiplicity Host;RecHit Multiplicity Device";
  layerCluster_RecHitMultiplicity_HostvsDevice_ = ibooker.book2I(histo, histoAxis, 100, 0, 100, 100, 0, 100);

  histo = "layerCluster_Layer_HostvsDevice";
  histoAxis = "layerCluster_Layer_HostvsDevice;Cluster Layer Host;Cluster Layer Device";
  layerCluster_Layer_HostvsDevice_ = ibooker.book2I(histo, histoAxis, 4, 0, 3, 4, 0, 3);

  histo = "layerCluster_Depth_HostvsDevice";
  histoAxis = "layerCluster_Depth_HostvsDevice;Cluster Depth Host;Cluster Depth Device";
  layerCluster_Depth_HostvsDevice_ = ibooker.book2I(histo, histoAxis, 8, 0, 7, 8, 0, 7);

  histo = "layerCluster_Eta_HostvsDevice";
  histoAxis = "layerCluster_Eta_HostvsDevice;Cluster #eta Host;Cluster #eta Device";
  layerCluster_Eta_HostvsDevice_ = ibooker.book2D(histo, histoAxis, 100, -5.f, 5.f, 100, -5.f, 5.f);

  histo = "layerCluster_Phi_HostvsDevice";
  histoAxis = "layerCluster_Phi_HostvsDevice;Cluster #phi Host;Cluster #phi Device";
  layerCluster_Phi_HostvsDevice_ = ibooker.book2D(histo, histoAxis, 100, -M_PI, M_PI, 100, -M_PI, M_PI);

  histo = "layerCluster_DuplicateMatches_HostvsDevice";
  histoAxis = "layerCluster_Duplicates_HostvsDevice;Cluster Duplicates Host;Cluster Duplicates Device";
  layerCluster_DuplicateMatches_HostvsDevice_ = ibooker.book1I(histo, histoAxis, 100, 0., 1000);

  layerCluster_Multiplicity_Diff_HostvsDevice_ = ibooker.book1D(
      "MultiplicityDiff", "PFCluster Multiplicity Difference; (Reference - Target);#entries", 100, -2, 2);
  layerCluster_Energy_Diff_HostvsDevice_ =
      ibooker.book1D("EnergyDiff", "PFCluster Energy Difference; (Reference - Target);#entries", 100, -2, 2);
  layerCluster_RecHitMultiplicity_Diff_HostvsDevice_ = ibooker.book1D(
      "RHMultiplicityDiff", "PFCluster RecHit Multiplicity Difference; (Reference - Target);#entries", 100, -2, 2);
  layerCluster_Layer_Diff_HostvsDevice_ =
      ibooker.book1D("LayerDiff", "PFCluster Layer Difference; (Reference - Target);#entries", 100, -2, 2);
  layerCluster_Depth_Diff_HostvsDevice_ =
      ibooker.book1D("DepthDiff", "PFCluster Depth Difference; (Reference - Target);#entries", 100, -2, 2);
  layerCluster_Eta_Diff_HostvsDevice_ =
      ibooker.book1D("EtaDiff", "PFCluster #eta Difference; (Reference - Target);#entries", 100, -0.5, 0.5);
  layerCluster_Phi_Diff_HostvsDevice_ =
      ibooker.book1D("PhiDiff", "PFCluster #phi Difference; (Reference - Target);#entries", 100, -0.5, 0.5);
}

void HGCalLayerClustersGPUComparisonTask::_resetMonitors(hcaldqm::UpdateFreq uf) { DQTask::_resetMonitors(uf); }

void HGCalLayerClustersGPUComparisonTask::_process(edm::Event const& event, edm::EventSetup const&) {
  const auto& layerClusters_ref = event.getHandle(layerClusterTok_ref_);
  const auto& layerClusters_target = event.getHandle(layerClusterTok_target_);

  // Exit early if any handle is invalid
  if (!layerClusters_ref || !layerClusters_target) {
    edm::LogWarning out("HGCalLayerClustersGPUComparisonTask");
    if (!layerClusters_ref)
      out << "reference PF cluster collection not found; ";
    if (!layerClusters_target)
      out << "target PF cluster collection not found; ";
    out << "the comparison will not run.";
    return;
  }

  auto lumiCache = luminosityBlockCache(event.getLuminosityBlock().index());
  _currentLS = lumiCache->currentLS;
  // Compare per-event PF cluster multiplicity

  if (layerClusters_ref->size() != layerClusters_target->size())
    LOGVERB("PFCaloGPUComparisonTask") << " PFCluster multiplicity " << layerClusters_ref->size() << " "
                                       << layerClusters_target->size();
  layerCluster_Multiplicity_HostvsDevice_->Fill((float)layerClusters_ref->size(), (float)layerClusters_ref->size());
  layerCluster_Multiplicity_Diff_HostvsDevice_->Fill((float)layerClusters_ref->size() - (float)layerClusters_target->size());
  //
  // Find matching PF cluster pairs
  std::vector<int> matched_idx;
  matched_idx.reserve(layerClusters_ref->size());
  for (unsigned i = 0; i < layerClusters_ref->size(); ++i) {
    bool matched = false;
    for (unsigned j = 0; j < layerClusters_target->size(); ++j) {
      if (layerClusters_ref->at(i).seed() == layerClusters_target->at(j).seed()) {
        if (!matched) {
          matched = true;
          matched_idx.push_back((int)j);
        } else {
          edm::LogWarning("PFCaloGPUComparisonTask") << "Found duplicate match";
          layerCluster_DuplicateMatches_HostvsDevice_->Fill((int)j);
        }
      }
    }
    if (!matched)
      matched_idx.push_back(-1);  // if you don't find a match, put a dummy number
  }

  //
  // Plot matching PF cluster variables
  for (unsigned i = 0; i < layerClusters_ref->size(); ++i) {
    if (matched_idx[i] >= 0) {
      unsigned int j = matched_idx[i];
      int ref_energy_bin =
          layerCluster_Energy_HostvsDevice_->getTH2F()->GetXaxis()->FindBin(layerClusters_ref->at(i).energy());
      int target_energy_bin =
          layerCluster_Energy_HostvsDevice_->getTH2F()->GetXaxis()->FindBin(layerClusters_target->at(j).energy());
      if (ref_energy_bin != target_energy_bin)
        edm::LogPrint("PFCaloGPUComparisonTask")
            << "Off-diagonal energy bin entries: " << layerClusters_ref->at(i).energy() << " "
            << layerClusters_ref->at(i).eta() << " " << layerClusters_ref->at(i).phi() << " "
            << layerClusters_target->at(j).energy() << " " << layerClusters_target->at(j).eta() << " "
            << layerClusters_target->at(j).phi() << std::endl;
      layerCluster_Energy_HostvsDevice_->Fill(layerClusters_ref->at(i).energy(), layerClusters_target->at(j).energy());
      layerCluster_Eta_HostvsDevice_->Fill(layerClusters_ref->at(i).eta(), layerClusters_target->at(j).eta());
      layerCluster_Phi_HostvsDevice_->Fill(layerClusters_ref->at(i).phi(), layerClusters_target->at(j).phi());
      layerCluster_Energy_Diff_HostvsDevice_->Fill(layerClusters_ref->at(i).energy() - layerClusters_target->at(j).energy());
      layerCluster_Eta_Diff_HostvsDevice_->Fill(layerClusters_ref->at(i).eta() - layerClusters_target->at(j).eta());
      layerCluster_Phi_Diff_HostvsDevice_->Fill(
          reco::deltaPhi(layerClusters_ref->at(i).phi(), layerClusters_target->at(j).phi()));
    }
  }
}

std::shared_ptr<hcaldqm::Cache> HGCalLayerClustersGPUComparisonTask::globalBeginLuminosityBlock(edm::LuminosityBlock const& lb,
                                                                                    edm::EventSetup const& es) const {
  return DQTask::globalBeginLuminosityBlock(lb, es);
}

void HGCalLayerClustersGPUComparisonTask::globalEndLuminosityBlock(edm::LuminosityBlock const& lb, edm::EventSetup const& es) {
  if (_ptype != fOnline)
    return;

  auto lumiCache = luminosityBlockCache(lb.index());
  _currentLS = lumiCache->currentLS;

  //	in the end always do the DQTask::endLumi
  DQTask::globalEndLuminosityBlock(lb, es);
}

void HGCalLayerClustersGPUComparisonTask::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<std::string>("subsystem", "ParticleFlow");
  desc.addUntracked<std::string>("name", "ParticleFlow/pfCaloGPUCompDir");
  desc.addUntracked<edm::InputTag>("layerClusterToken_ref", edm::InputTag("hltParticleFlowClusterHCALSerialSync"));
  desc.addUntracked<edm::InputTag>("layerClusterToken_target", edm::InputTag("hltParticleFlowClusterHCAL"));
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(HGCalLayerClustersGPUComparisonTask);
