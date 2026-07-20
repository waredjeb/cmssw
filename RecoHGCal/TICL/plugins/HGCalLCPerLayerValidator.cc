// HGCalLCPerLayerValidator
// ------------------------------------------------------------------------
// Lightweight per-agent validator for the PataTune LayerClusters optimization.
// For one (cloned) merged LayerCluster collection + its LC<->CaloParticle
// association, it accumulates the counters needed to build three objectives:
//
//   * per-(CaloParticle, layer) efficiency
//   * number of LayerClusters
//   * per-(CaloParticle, layer) splitting (over-segmentation)
//
// One TTree entry is filled per event; the driver sums the columns with uproot.
// TFileService namespaces the tree under this module's label, so each agent's
// clone writes into its own directory (e.g. lcValidatorAgent0/output).
//
// Association score conventions (matching HGCalValidator):
//   RecoToSim quality = float score;              LC "pure" if score < scoreCutLCtoCP
//   SimToReco quality = pair<float,float>, .second = score;  CP<-LC if .second < scoreCutCPtoLC

#include <set>
#include <map>
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/Associations/interface/LayerClusterToCaloParticleAssociator.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "TTree.h"

class HGCalLCPerLayerValidator : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit HGCalLCPerLayerValidator(const edm::ParameterSet&);
  ~HGCalLCPerLayerValidator() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  static bool isHGCal(const DetId& id) {
    const int d = id.det();
    return d == DetId::HGCalEE || d == DetId::HGCalHSi || d == DetId::HGCalHSc || d == DetId::Forward;
  }

  const edm::EDGetTokenT<reco::CaloClusterCollection> lcToken_;
  const edm::EDGetTokenT<std::vector<CaloParticle>> cpToken_;
  const edm::EDGetTokenT<ticl::RecoToSimCollectionT<reco::CaloClusterCollection>> recoToSimToken_;
  const edm::EDGetTokenT<ticl::SimToRecoCollectionT<reco::CaloClusterCollection>> simToRecoToken_;
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;

  const double scoreCutLCtoCP_;
  const double scoreCutCPtoLC_;

  std::shared_ptr<hgcal::RecHitTools> tools_;

  // per-event branch buffers (all summed by the driver)
  TTree* tree_;
  unsigned int b_nCP_;
  unsigned int b_nLC_;
  unsigned int b_nCPLayer_total_;   // (CP,layer) pairs with sim energy   -> eff denominator
  unsigned int b_nCPLayer_eff_;     // ... with >=1 associated LC on layer -> eff numerator
  unsigned int b_nCPLayer_split_;   // sum max(0, #assocLC_on_layer - 1)
  unsigned int b_nLC_merged_;       // LCs associated to >1 CP (score<cut)
};

HGCalLCPerLayerValidator::HGCalLCPerLayerValidator(const edm::ParameterSet& ps)
    : lcToken_(consumes<reco::CaloClusterCollection>(ps.getParameter<edm::InputTag>("layerClusters"))),
      cpToken_(consumes<std::vector<CaloParticle>>(ps.getParameter<edm::InputTag>("caloParticles"))),
      recoToSimToken_(consumes<ticl::RecoToSimCollectionT<reco::CaloClusterCollection>>(
          ps.getParameter<edm::InputTag>("lcToCpAssociator"))),
      simToRecoToken_(consumes<ticl::SimToRecoCollectionT<reco::CaloClusterCollection>>(
          ps.getParameter<edm::InputTag>("lcToCpAssociator"))),
      caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
      scoreCutLCtoCP_(ps.getParameter<double>("scoreCutLCtoCP")),
      scoreCutCPtoLC_(ps.getParameter<double>("scoreCutCPtoLC")) {
  usesResource(TFileService::kSharedResource);
  tools_ = std::make_shared<hgcal::RecHitTools>();

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("output", "per-event LC tuning counters");
  tree_->Branch("nCP", &b_nCP_);
  tree_->Branch("nLC", &b_nLC_);
  tree_->Branch("nCPLayer_total", &b_nCPLayer_total_);
  tree_->Branch("nCPLayer_eff", &b_nCPLayer_eff_);
  tree_->Branch("nCPLayer_split", &b_nCPLayer_split_);
  tree_->Branch("nLC_merged", &b_nLC_merged_);
}

void HGCalLCPerLayerValidator::analyze(const edm::Event& evt, const edm::EventSetup& es) {
  tools_->setGeometry(es.getData(caloGeomToken_));

  const auto& lcs = evt.get(lcToken_);
  const auto lcHandle = evt.getHandle(lcToken_);
  const auto cpHandle = evt.getHandle(cpToken_);
  const auto& recoToSim = evt.get(recoToSimToken_);
  const auto& simToReco = evt.get(simToRecoToken_);

  b_nCP_ = 0;
  b_nLC_ = lcs.size();
  b_nCPLayer_total_ = 0;
  b_nCPLayer_eff_ = 0;
  b_nCPLayer_split_ = 0;
  b_nLC_merged_ = 0;

  // ---- LC layer lookup (offset layer: EE 1..N, then FH, then BH) ----
  std::vector<unsigned int> lcLayer(lcs.size(), 0u);
  for (size_t i = 0; i < lcs.size(); ++i) {
    DetId seed = lcs[i].seed();
    if (seed.rawId() == 0 && !lcs[i].hitsAndFractions().empty())
      seed = lcs[i].hitsAndFractions().front().first;
    if (seed.rawId() != 0 && isHGCal(seed))
      lcLayer[i] = tools_->getLayerWithOffset(seed);
  }

  // ---- merge: LCs associated (score<cut) to more than one CaloParticle ----
  for (size_t i = 0; i < lcs.size(); ++i) {
    edm::Ref<reco::CaloClusterCollection> lcRef(lcHandle, i);
    auto it = recoToSim.find(lcRef);
    if (it == recoToSim.end())
      continue;
    unsigned int nCpMatched = 0;
    for (const auto& cpWithScore : it->val) {
      if (cpWithScore.second < scoreCutLCtoCP_)
        ++nCpMatched;
    }
    if (nCpMatched > 1)
      ++b_nLC_merged_;
  }

  // ---- per-(CP, layer) efficiency + splitting ----
  const auto& caloParticles = *cpHandle;
  for (size_t c = 0; c < caloParticles.size(); ++c) {
    const CaloParticle& cp = caloParticles[c];

    // layers where this CaloParticle deposits sim energy
    std::set<unsigned int> cpLayers;
    for (const auto& hf : cp.hits_and_fractions()) {
      DetId id(hf.first);
      if (isHGCal(id))
        cpLayers.insert(tools_->getLayerWithOffset(id));
    }
    if (cpLayers.empty())
      continue;
    ++b_nCP_;
    b_nCPLayer_total_ += cpLayers.size();

    // count well-associated LCs per layer for this CP
    std::map<unsigned int, unsigned int> assocLCByLayer;
    CaloParticleRef cpRef(cpHandle, c);
    auto it = simToReco.find(cpRef);
    if (it != simToReco.end()) {
      for (const auto& lcWithQual : it->val) {
        const double score = lcWithQual.second.second;  // pair<sharedE, score>
        if (score < scoreCutCPtoLC_)
          ++assocLCByLayer[lcLayer[lcWithQual.first.index()]];
      }
    }

    for (unsigned int layer : cpLayers) {
      auto ait = assocLCByLayer.find(layer);
      const unsigned int n = (ait == assocLCByLayer.end()) ? 0u : ait->second;
      if (n >= 1)
        ++b_nCPLayer_eff_;
      if (n > 1)
        b_nCPLayer_split_ += (n - 1);
    }
  }

  tree_->Fill();
}

void HGCalLCPerLayerValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("caloParticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("lcToCpAssociator", edm::InputTag("layerClusterCaloParticleAssociation"));
  desc.add<double>("scoreCutLCtoCP", 0.1);
  desc.add<double>("scoreCutCPtoLC", 0.1);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(HGCalLCPerLayerValidator);
