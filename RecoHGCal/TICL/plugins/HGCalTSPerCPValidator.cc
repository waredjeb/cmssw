// HGCalTSPerCPValidator
// ------------------------------------------------------------------------
// Lightweight per-agent trackster validator for the PataTune joint
// CLUE2D+CLUE3D optimization. For one (cloned) trackster collection, its
// (cloned) merged LayerClusters and the corresponding LC<->CaloParticle
// association, it accumulates the counters needed to build three objectives:
//
//   * efficiency  : CaloParticles with >=1 trackster at shared energy
//                   fraction >= minSharedEneFracEfficiency (0.5)
//   * split       : CaloParticles shared by >1 trackster at shared energy
//                   fraction >= minSharedEneFracSplit (0.1) -- over-segmentation
//   * merge       : tracksters matching >1 CaloParticle at scoreCutRecoToSim
//                   (over-clustering)
//
// Why split is defined on shared energy and not on the score
// ----------------------------------------------------------
// The official "duplicate" counter (sim->reco score < 0.2) is blind to
// splitting, and so is any score-based variant. With simFraction ~ 1 the
// sim->reco score of a trackster covering a subset S of the CaloParticle's
// layer clusters reduces to
//     score ~ 1 - sum_{lc in S} E^2 / sum_{all lc} E^2
// i.e. one minus the squared-energy fraction that piece carries. A CP broken
// in two therefore gives ~0.5 per piece, a CP broken in three gives ~0.67, and
// so on: raising the cut to catch 2-way splits still misses 3-way splits, and
// the counter is NON-MONOTONIC in the severity of the splitting -- worse
// fragmentation reads as less splitting once every piece falls past the cut.
// Measured on a 4-particle CloseBy sample: ~3.5 tracksters per CaloParticle
// with zero fakes, yet the score-based duplicate counter was identically 0.
//
// Counting pieces by shared energy fraction instead is monotonic (an n-way
// split yields n-1 excess pieces) and reuses exactly the quantity the
// efficiency is already built on: efficiency asks "is there one piece
// >= 0.5?", split asks "how many pieces >= 0.1?".
//
// The score-based purity/duplicate counters are kept as diagnostic branches.
//
// One TTree entry is filled per event; the driver sums the columns with uproot.
// TFileService namespaces the tree under this module's label, so each agent's
// clone writes into its own directory (e.g. tsValidatorAgent0/output). This
// mirrors HGCalLCPerLayerValidator exactly.
//
// Why no SimTracksters
// --------------------
// The official trackster<->SimTrackster association would require
// ticlSimTracksters, which drags in the whole tracking-truth chain
// (generalTracks, TrackingParticles, simHitTPAssocProducer,
// trackingParticleRecoTrackAsssociation) -- and, because the tuned CLUE 2D
// parameters move the LayerClusters, it would have to be cloned per agent.
//
// It is also unnecessary. SimTrackstersProducer::addTrackster builds the
// "fromCPs" sim tracksters purely from the LC<->CaloParticle association:
//     fraction            = sharedEnergy / lc.energy()
//     vertices           <- the associated LCs
//     vertex_multiplicity = 1 / fraction
// and AllLayerClusterToTracksterAssociatorsProducer then recovers
//     simSharedEnergy = lc.energy() / vertex_multiplicity = sharedEnergy
// i.e. exactly the shared energy already stored in the LC<->CP association.
// Tracks only decorate the sim trackster with trackIdx / regressed energy,
// which no score here uses. So we consume the LC<->CP association directly and
// apply the official score formulas verbatim.
//
// Score formulas are copied from
//   SimCalorimetry/HGCalAssociatorProducers/plugins/
//     AllTracksterToSimTracksterAssociatorsByLCsProducer.cc
// and the cut values from Validation/HGCalValidation/src/HGVHistoProducerAlgo.cc
//   ScoreCutTStoSTSFakeMerge_ = 0.6   (reco->sim, fake & merge)
//   ScoreCutSTStoTSPurDup_    = 0.2   (sim->reco, purity & duplicate)
// and HGCalValidator.cc  minTSTSharedEneFracEfficiency_ = 0.5.
//
// Deviations from the official numbers, deliberate and documented:
//   * no CaloParticle preselection. HGCalValidator.cc:349 applies
//     CaloParticleSelector (CaloParticleSelectionForEfficiency_cfi: ptMinCP=0.5,
//     ptMaxCP=300, |y|<3.1, notConvertedOnlyCP=True, pdgId whitelist); we count
//     every CP with a non-empty sim trackster, matching HGCalLCPerLayerValidator.
//   * no filteredLayerClustersSimTracksters input mask (that filter uses
//     min_cluster_size=0, so it accepts essentially every LC anyway).
//
// Validation gate, measured 100 events of the mixed gamma+pion CloseBy sample,
// this validator run in the SAME process as hgcalValidator on the SAME
// ticlTrackstersCLUE3DHigh collection (optimization/valgate/step4_val_grafted.py),
// compared against HGCalValidator/ticlTrackstersCLUE3DHigh/TSbyLCs_CP:
//
//                        official     ours
//   trackster denom          3129     3129   identical
//   merged tracksters         258      258   identical
//   merge rate             0.0825   0.0825   identical
//   fake rate              0.0000   0.0000   identical
//   CP denom                  387      398   +2.8%  (CaloParticleSelector)
//   efficiency             0.7468   0.7362   -1.1pp (same cause)
//   purity                 0.7416   0.7312   -1.0pp (same cause)
//   duplicate              0.0000   0.0000   identical -- dead on both sides
//
// The reco->sim side is exact. The sim->reco side differs only through the
// CaloParticle preselection above. Note the official duplicate counter is zero
// on the official chain too, which is what motivated the shared-energy split
// definition; the same events give split = 0.2236.

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

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

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/Associations/interface/LayerClusterToCaloParticleAssociator.h"

#include "TTree.h"

class HGCalTSPerCPValidator : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit HGCalTSPerCPValidator(const edm::ParameterSet&);
  ~HGCalTSPerCPValidator() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<std::vector<ticl::Trackster>> tsToken_;
  const edm::EDGetTokenT<reco::CaloClusterCollection> lcToken_;
  const edm::EDGetTokenT<std::vector<CaloParticle>> cpToken_;
  const edm::EDGetTokenT<ticl::SimToRecoCollectionT<reco::CaloClusterCollection>> simToRecoToken_;

  const double scoreCutRecoToSim_;  // fake / merge     (official 0.6)
  const double scoreCutSimToReco_;  // purity / duplicate (official 0.2, diagnostics)
  const double minSharedEneFracEfficiency_;  // official 0.5
  const double minSharedEneFracSplit_;       // 0.1

  // per-event branch buffers (all summed by the driver)
  TTree* tree_;
  unsigned int b_nCP_;          // CPs with a non-empty sim trackster -> eff/split denominator
  unsigned int b_nCP_eff_;      // ... with >=1 TS at sharedEneFrac >= effCut -> eff numerator
  unsigned int b_nCP_split_;    // ... with >1  TS at sharedEneFrac >= splitCut -> split numerator
  unsigned int b_nCP_split_excess_;  // sum max(0, #such TS - 1) -> splitting severity (diag)
  unsigned int b_nCP_pure_;     // ... with >=1 TS at simToReco score < cut (diag)
  unsigned int b_nCP_dup_;      // sum max(0, #TS passing score cut - 1)    (diag)
  unsigned int b_nTS_;          // tracksters with >=1 vertex -> merge denominator
  unsigned int b_nTS_fake_;     // ... matching ZERO CPs at recoToSim score < cut (diag)
  unsigned int b_nTS_merged_;   // ... matching >1 CP at recoToSim score < cut -> merge numerator
};

HGCalTSPerCPValidator::HGCalTSPerCPValidator(const edm::ParameterSet& ps)
    : tsToken_(consumes<std::vector<ticl::Trackster>>(ps.getParameter<edm::InputTag>("tracksters"))),
      lcToken_(consumes<reco::CaloClusterCollection>(ps.getParameter<edm::InputTag>("layerClusters"))),
      cpToken_(consumes<std::vector<CaloParticle>>(ps.getParameter<edm::InputTag>("caloParticles"))),
      simToRecoToken_(consumes<ticl::SimToRecoCollectionT<reco::CaloClusterCollection>>(
          ps.getParameter<edm::InputTag>("lcToCpAssociator"))),
      scoreCutRecoToSim_(ps.getParameter<double>("scoreCutRecoToSim")),
      scoreCutSimToReco_(ps.getParameter<double>("scoreCutSimToReco")),
      minSharedEneFracEfficiency_(ps.getParameter<double>("minSharedEneFracEfficiency")),
      minSharedEneFracSplit_(ps.getParameter<double>("minSharedEneFracSplit")) {
  usesResource(TFileService::kSharedResource);

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("output", "per-event trackster tuning counters");
  tree_->Branch("nCP", &b_nCP_);
  tree_->Branch("nCP_eff", &b_nCP_eff_);
  tree_->Branch("nCP_split", &b_nCP_split_);
  tree_->Branch("nCP_split_excess", &b_nCP_split_excess_);
  tree_->Branch("nCP_pure", &b_nCP_pure_);
  tree_->Branch("nCP_dup", &b_nCP_dup_);
  tree_->Branch("nTS", &b_nTS_);
  tree_->Branch("nTS_fake", &b_nTS_fake_);
  tree_->Branch("nTS_merged", &b_nTS_merged_);
}

void HGCalTSPerCPValidator::analyze(const edm::Event& evt, const edm::EventSetup&) {
  const auto& tracksters = evt.get(tsToken_);
  const auto& lcs = evt.get(lcToken_);
  const auto cpHandle = evt.getHandle(cpToken_);
  const auto& simToReco = evt.get(simToRecoToken_);
  const auto& caloParticles = *cpHandle;

  b_nCP_ = 0;
  b_nCP_eff_ = 0;
  b_nCP_split_ = 0;
  b_nCP_split_excess_ = 0;
  b_nCP_pure_ = 0;
  b_nCP_dup_ = 0;
  b_nTS_ = 0;
  b_nTS_fake_ = 0;
  b_nTS_merged_ = 0;

  // ------------------------------------------------------------------ sim ---
  // The "sim trackster" of CaloParticle c is the set of LayerClusters the
  // LC<->CP association attaches to it, with simSharedEnergy = the association's
  // shared energy (see the header comment for why this equals the official one).
  //
  // simVerts[c] : (lcIndex, simSharedEnergy) for CaloParticle c
  // lcToCp[lc]  : (cpSlot,  simSharedEnergy) for LayerCluster lc  (cpSlot indexes simVerts)
  std::vector<std::vector<std::pair<unsigned int, float>>> simVerts;
  std::vector<float> simRawEnergy;
  std::vector<std::vector<std::pair<unsigned int, float>>> lcToCp(lcs.size());

  simVerts.reserve(caloParticles.size());
  simRawEnergy.reserve(caloParticles.size());
  for (size_t c = 0; c < caloParticles.size(); ++c) {
    CaloParticleRef cpRef(cpHandle, c);
    auto it = simToReco.find(cpRef);
    if (it == simToReco.end())
      continue;

    std::vector<std::pair<unsigned int, float>> verts;
    float rawEnergy = 0.f;
    for (const auto& lcWithQual : it->val) {
      const unsigned int lcIdx = lcWithQual.first.index();
      const float sharedEnergy = lcWithQual.second.first;  // pair<sharedE, score>
      if (lcIdx >= lcs.size() || lcs[lcIdx].energy() <= 0.f || sharedEnergy <= 0.f)
        continue;
      verts.emplace_back(lcIdx, sharedEnergy);
      rawEnergy += sharedEnergy;
    }
    if (verts.empty() || rawEnergy <= 0.f)
      continue;

    const unsigned int cpSlot = simVerts.size();
    for (const auto& [lcIdx, sharedEnergy] : verts)
      lcToCp[lcIdx].emplace_back(cpSlot, sharedEnergy);
    simVerts.push_back(std::move(verts));
    simRawEnergy.push_back(rawEnergy);
  }
  b_nCP_ = simVerts.size();

  // ----------------------------------------------------------------- reco ---
  // recoVerts[t] : (lcIndex, recoFraction) for trackster t
  // lcToTs[lc]   : (tsSlot,  recoSharedEnergy) for LayerCluster lc
  std::vector<std::vector<std::pair<unsigned int, float>>> recoVerts;
  std::vector<std::vector<std::pair<unsigned int, float>>> lcToTs(lcs.size());

  recoVerts.reserve(tracksters.size());
  for (const auto& ts : tracksters) {
    std::vector<std::pair<unsigned int, float>> verts;
    for (unsigned int i = 0; i < ts.vertices().size(); ++i) {
      const unsigned int lcIdx = ts.vertices()[i];
      const float mult = ts.vertex_multiplicity(i);
      if (lcIdx >= lcs.size() || lcs[lcIdx].energy() <= 0.f || mult <= 0.f)
        continue;
      verts.emplace_back(lcIdx, 1.f / mult);
    }
    if (verts.empty())
      continue;

    const unsigned int tsSlot = recoVerts.size();
    for (const auto& [lcIdx, recoFraction] : verts)
      lcToTs[lcIdx].emplace_back(tsSlot, lcs[lcIdx].energy() * recoFraction);
    recoVerts.push_back(std::move(verts));
  }
  b_nTS_ = recoVerts.size();

  // ------------------------------------------- reco -> sim : fake & merge ---
  for (const auto& verts : recoVerts) {
    float denominator = 0.f;
    std::vector<unsigned int> assocCps;
    for (const auto& [lcIdx, recoFraction] : verts) {
      const float e = lcs[lcIdx].energy();
      denominator += e * e * recoFraction * recoFraction;
      for (const auto& [cpSlot, simSharedEnergy] : lcToCp[lcIdx])
        assocCps.push_back(cpSlot);
    }
    if (denominator <= 0.f)
      continue;

    std::sort(assocCps.begin(), assocCps.end());
    assocCps.erase(std::unique(assocCps.begin(), assocCps.end()), assocCps.end());
    if (assocCps.empty()) {
      ++b_nTS_fake_;
      continue;
    }

    const float invDenominator = 1.f / denominator;
    std::unordered_map<unsigned int, float> scoreByCp;
    for (const auto& [lcIdx, recoFraction] : verts) {
      const float e = lcs[lcIdx].energy();
      const float invE = 1.f / e;
      // CaloParticles associated to this trackster but NOT sharing this LC
      // contribute with simFraction = 0, i.e. a full recoFraction penalty.
      for (unsigned int cpSlot : assocCps) {
        float simSharedEnergy = 0.f;
        for (const auto& [slot, sharedEnergy] : lcToCp[lcIdx]) {
          if (slot == cpSlot) {
            simSharedEnergy = sharedEnergy;
            break;
          }
        }
        const float simFraction = simSharedEnergy * invE;
        const float d = std::max(0.f, recoFraction - simFraction);
        scoreByCp[cpSlot] += invDenominator * d * d * e * e;
      }
    }

    unsigned int nCpMatched = 0;
    for (const auto& [cpSlot, score] : scoreByCp) {
      if (score < scoreCutRecoToSim_)
        ++nCpMatched;
    }
    if (nCpMatched == 0)
      ++b_nTS_fake_;
    if (nCpMatched > 1)
      ++b_nTS_merged_;
  }

  // -------------------------------------- sim -> reco : efficiency & dup ---
  for (unsigned int cpSlot = 0; cpSlot < simVerts.size(); ++cpSlot) {
    const auto& verts = simVerts[cpSlot];

    float denominator = 0.f;
    std::vector<unsigned int> assocTs;
    for (const auto& [lcIdx, simSharedEnergy] : verts) {
      denominator += simSharedEnergy * simSharedEnergy;
      for (const auto& [tsSlot, recoSharedEnergy] : lcToTs[lcIdx])
        assocTs.push_back(tsSlot);
    }
    if (denominator <= 0.f)
      continue;

    std::sort(assocTs.begin(), assocTs.end());
    assocTs.erase(std::unique(assocTs.begin(), assocTs.end()), assocTs.end());
    if (assocTs.empty())
      continue;

    const float invDenominator = 1.f / denominator;
    std::unordered_map<unsigned int, float> scoreByTs;
    std::unordered_map<unsigned int, float> sharedByTs;
    for (const auto& [lcIdx, simSharedEnergy] : verts) {
      const float e = lcs[lcIdx].energy();
      const float invE = 1.f / e;
      const float simFraction = simSharedEnergy * invE;
      // Tracksters associated to this CaloParticle but NOT sharing this LC
      // contribute with recoFraction = 0, i.e. a full simFraction penalty.
      for (unsigned int tsSlot : assocTs) {
        float recoSharedEnergy = 0.f;
        for (const auto& [slot, sharedEnergy] : lcToTs[lcIdx]) {
          if (slot == tsSlot) {
            recoSharedEnergy = sharedEnergy;
            break;
          }
        }
        const float recoFraction = recoSharedEnergy * invE;
        const float d = std::max(0.f, simFraction - recoFraction);
        scoreByTs[tsSlot] += invDenominator * d * d * e * e;
        sharedByTs[tsSlot] += std::min(simSharedEnergy, recoSharedEnergy);
      }
    }

    const float invSimRawEnergy = 1.f / simRawEnergy[cpSlot];
    bool efficient = false;
    unsigned int nTsShared = 0;   // pieces carrying a real fraction of this CP
    unsigned int nTsMatched = 0;  // official score-based purity/duplicate (diagnostic)
    for (const auto& [tsSlot, score] : scoreByTs) {
      const float sharedFraction = sharedByTs[tsSlot] * invSimRawEnergy;
      if (sharedFraction >= minSharedEneFracEfficiency_)
        efficient = true;
      if (sharedFraction >= minSharedEneFracSplit_)
        ++nTsShared;
      if (score < scoreCutSimToReco_)
        ++nTsMatched;
    }
    if (efficient)
      ++b_nCP_eff_;
    if (nTsShared > 1) {
      ++b_nCP_split_;
      b_nCP_split_excess_ += (nTsShared - 1);
    }
    if (nTsMatched > 0)
      ++b_nCP_pure_;
    if (nTsMatched > 1)
      b_nCP_dup_ += (nTsMatched - 1);
  }

  tree_->Fill();
}

void HGCalTSPerCPValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracksters", edm::InputTag("ticlTrackstersCLUE3DHigh"));
  desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("caloParticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("lcToCpAssociator", edm::InputTag("layerClusterCaloParticleAssociation"));
  desc.add<double>("scoreCutRecoToSim", 0.6);
  desc.add<double>("scoreCutSimToReco", 0.2);
  desc.add<double>("minSharedEneFracEfficiency", 0.5);
  desc.add<double>("minSharedEneFracSplit", 0.1);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(HGCalTSPerCPValidator);
