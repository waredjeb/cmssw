// DQM histograms of the particle-flow candidate truth validation, filled from the two
// records PFCandidateTruthAssociator persists. No physics decision is taken here: every
// flag and fraction comes from the records, so the validator only counts.
//
// Folder layout, under <dirName>/<collection key>/:
//   truth/<species>/<region>/
//     num_simul_{pt,eta,caloeta}          every particle of the level
//     num_<rung>_expected_{pt,eta,caloeta} particles the rung applies to (its denominator)
//     num_<rung>_{pt,eta,caloeta}          of those, the ones that pass
//     ladder_trackFound, ladder_trackMissing   the cut-flow per branch of the tree: bin k
//                                          counts particles passing rungs 1..k of that
//                                          branch, bin 0 the branch population
//     response_trackFound, response_trackMissing   candidate over truth energy
//     foreign_fraction                     foreign energy share of the resolved candidate
//   reco/<pftype>/<region>/
//     num_reco_{pt,eta}, num_<class>_{pt,eta}, own_fraction
// <region> is one of barrel, transition, endcap, forward, all; <species> and <pftype>
// include an "all" entry as well.

#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "SimDataFormats/TruthInfo/interface/PFCandidateTruthRecords.h"

namespace {
  using truth::PFCandidateClass;
  using truth::PFRung;

  const std::vector<std::string> kRegions = {"barrel", "transition", "endcap", "forward"};
  const std::vector<std::string> kSpecies = {"pi", "K", "p", "e", "mu", "gamma", "pi0", "n", "K0L", "other"};
  const std::vector<std::string> kPFTypes = {"X", "h", "e", "mu", "gamma", "h0", "h_HF", "egamma_HF"};
  const std::vector<std::string> kClasses = {"unmatched", "matched", "merged", "split", "other", "unevaluated"};

  std::string speciesOf(int pdgId) {
    switch (std::abs(pdgId)) {
      case 211:
        return "pi";
      case 321:
        return "K";
      case 2212:
        return "p";
      case 11:
        return "e";
      case 13:
        return "mu";
      case 22:
        return "gamma";
      case 111:
        return "pi0";
      case 2112:
        return "n";
      case 130:
        return "K0L";
      default:
        return "other";
    }
  }

  // The rungs, each with the flag that says it applies (its denominator) and the flag
  // that says it passed. CandidateFound and the candidate rungs apply to every particle.
  struct RungDef {
    const char* name;
    PFRung expected;
    PFRung passed;
    bool alwaysExpected;
  };
  const std::vector<RungDef> kRungs = {
      {"trackFound", PFRung::TrackExpected, PFRung::TrackFound, false},
      {"trackInCandidate", PFRung::TrackFound, PFRung::TrackInCandidate, false},
      {"ecalCollected", PFRung::EcalExpected, PFRung::EcalCollected, false},
      {"ecalLinked", PFRung::EcalCollected, PFRung::EcalLinked, false},
      {"hcalCollected", PFRung::HcalExpected, PFRung::HcalCollected, false},
      {"hcalLinked", PFRung::HcalCollected, PFRung::HcalLinked, false},
      {"hgcalCollected", PFRung::HgcalExpected, PFRung::HgcalCollected, false},
      {"hgcalLinked", PFRung::HgcalCollected, PFRung::HgcalLinked, false},
      {"caloSameCandidate", PFRung::HcalCollected, PFRung::CaloSameCandidate, false},  // and ECAL collected, see fill
      {"candidateFound", PFRung::TrackExpected, PFRung::CandidateFound, true},
      {"merged", PFRung::CandidateFound, PFRung::Merged, false},
      {"partial", PFRung::CandidateFound, PFRung::Partial, false},
      {"clean", PFRung::CandidateFound, PFRung::Clean, false},
      {"pdgCorrect", PFRung::PdgEvaluated, PFRung::PdgCorrect, false},
  };

  // The cut-flow of each branch of the tree: rung k applies only if every previous one
  // passed (or was not expected, in which case it is skipped, not failed).
  const std::vector<PFRung> kLadderTrackFound = {PFRung::TrackFound,
                                                 PFRung::TrackInCandidate,
                                                 PFRung::EcalCollected,
                                                 PFRung::EcalLinked,
                                                 PFRung::HcalCollected,
                                                 PFRung::HcalLinked,
                                                 PFRung::HgcalCollected,
                                                 PFRung::HgcalLinked,
                                                 PFRung::CaloSameCandidate,
                                                 PFRung::Clean,
                                                 PFRung::PdgCorrect};
  const std::vector<const char*> kLadderTrackFoundLabels = {"population",
                                                            "track found",
                                                            "track in cand.",
                                                            "ECAL collected",
                                                            "ECAL linked",
                                                            "HCAL collected",
                                                            "HCAL linked",
                                                            "HGCAL collected",
                                                            "HGCAL linked",
                                                            "calo same cand.",
                                                            "clean",
                                                            "pdg correct"};
  const std::vector<PFRung> kLadderTrackMissing = {PFRung::EcalCollected,
                                                   PFRung::HcalCollected,
                                                   PFRung::HgcalCollected,
                                                   PFRung::CaloSameCandidate,
                                                   PFRung::CandidateFound,
                                                   PFRung::Clean,
                                                   PFRung::PdgCorrect};
  const std::vector<const char*> kLadderTrackMissingLabels = {"population",
                                                              "ECAL collected",
                                                              "HCAL collected",
                                                              "HGCAL collected",
                                                              "calo same cand.",
                                                              "candidate found",
                                                              "clean",
                                                              "pdg correct"};

  // Whether a rung of a ladder applies to this particle, from its Expected flags.
  bool rungApplies(truth::PFCandidateTruthRecord const& r, PFRung rung) {
    switch (rung) {
      case PFRung::TrackFound:
      case PFRung::TrackInCandidate:
        return r.has(PFRung::TrackExpected);
      case PFRung::EcalCollected:
      case PFRung::EcalLinked:
        return r.has(PFRung::EcalExpected);
      case PFRung::HcalCollected:
      case PFRung::HcalLinked:
        return r.has(PFRung::HcalExpected);
      case PFRung::HgcalCollected:
      case PFRung::HgcalLinked:
        return r.has(PFRung::HgcalExpected);
      case PFRung::CaloSameCandidate:
        return r.has(PFRung::EcalCollected) && r.has(PFRung::HcalCollected);
      case PFRung::PdgCorrect:
        return r.has(PFRung::PdgEvaluated);
      default:
        return true;
    }
  }
}  // namespace

class PFCandidateTruthValidator : public DQMEDAnalyzer {
public:
  explicit PFCandidateTruthValidator(edm::ParameterSet const&);
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  void analyze(edm::Event const&, edm::EventSetup const&) override;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  struct Axes {
    MonitorElement* pt = nullptr;
    MonitorElement* eta = nullptr;
    MonitorElement* caloeta = nullptr;
    void fill(double ptValue, double etaValue, double caloEtaValue) const {
      pt->Fill(ptValue);
      eta->Fill(etaValue);
      if (caloeta != nullptr)
        caloeta->Fill(caloEtaValue);
    }
  };
  struct TruthFolder {
    Axes simul;
    std::vector<Axes> expected;  // per rung
    std::vector<Axes> passed;    // per rung
    MonitorElement* ladderTrackFound = nullptr;
    MonitorElement* ladderTrackMissing = nullptr;
    MonitorElement* responseTrackFound = nullptr;
    MonitorElement* responseTrackMissing = nullptr;
    MonitorElement* foreignFraction = nullptr;
  };
  struct RecoFolder {
    Axes reco;
    std::vector<Axes> classes;
    MonitorElement* ownFraction = nullptr;
  };

  Axes bookAxes(DQMStore::IBooker& booker, std::string const& name, bool withCaloEta) const;
  void fillTruth(TruthFolder& folder, truth::PFCandidateTruthRecord const& r) const;
  void fillReco(RecoFolder& folder, truth::PFCandidateRecord const& c) const;

  const edm::EDGetTokenT<truth::PFCandidateTruthRecordCollection> truthToken_;
  const edm::EDGetTokenT<truth::PFCandidateRecordCollection> recoToken_;
  const std::string dirName_;
  std::vector<float> ptEdges_;
  const int nEta_;
  const double etaMax_;

  std::map<std::string, std::map<std::string, TruthFolder>> truth_;  // species -> region
  std::map<std::string, std::map<std::string, RecoFolder>> reco_;    // pftype -> region
};

PFCandidateTruthValidator::PFCandidateTruthValidator(edm::ParameterSet const& cfg)
    : truthToken_(consumes<truth::PFCandidateTruthRecordCollection>(cfg.getParameter<edm::InputTag>("truthRecords"))),
      recoToken_(consumes<truth::PFCandidateRecordCollection>(cfg.getParameter<edm::InputTag>("candidateRecords"))),
      dirName_(cfg.getParameter<std::string>("dirName")),
      nEta_(cfg.getParameter<int>("nintEta")),
      etaMax_(cfg.getParameter<double>("maxEta")) {
  // Logarithmic pt axis: a flat-pt gun and a ttbar spectrum are both readable on it.
  const int nPt = cfg.getParameter<int>("nintPt");
  const double lo = std::log10(cfg.getParameter<double>("minPt"));
  const double hi = std::log10(cfg.getParameter<double>("maxPt"));
  for (int i = 0; i <= nPt; ++i)
    ptEdges_.push_back(static_cast<float>(std::pow(10., lo + (hi - lo) * i / nPt)));
}

PFCandidateTruthValidator::Axes PFCandidateTruthValidator::bookAxes(DQMStore::IBooker& booker,
                                                                    std::string const& name,
                                                                    bool withCaloEta) const {
  Axes axes;
  axes.pt = booker.book1D(name + "_pt", name + ";p_{T} [GeV];", static_cast<int>(ptEdges_.size()) - 1, ptEdges_.data());
  axes.eta = booker.book1D(name + "_eta", name + ";#eta;", nEta_, -etaMax_, etaMax_);
  if (withCaloEta)
    axes.caloeta = booker.book1D(name + "_caloeta", name + ";#eta at calorimeter entrance;", nEta_, -etaMax_, etaMax_);
  return axes;
}

void PFCandidateTruthValidator::bookHistograms(DQMStore::IBooker& booker, edm::Run const&, edm::EventSetup const&) {
  std::vector<std::string> regions = kRegions;
  regions.push_back("all");
  std::vector<std::string> species = kSpecies;
  species.push_back("all");
  for (auto const& s : species) {
    for (auto const& region : regions) {
      booker.setCurrentFolder(dirName_ + "/truth/" + s + "/" + region);
      TruthFolder& f = truth_[s][region];
      f.simul = bookAxes(booker, "num_simul", true);
      for (auto const& rung : kRungs) {
        f.expected.push_back(bookAxes(booker, std::string("num_") + rung.name + "_expected", true));
        f.passed.push_back(bookAxes(booker, std::string("num_") + rung.name, true));
      }
      auto bookLadder = [&](const char* name, std::vector<const char*> const& labels) {
        MonitorElement* me = booker.book1D(name, std::string(name) + ";;particles", labels.size(), 0., labels.size());
        for (std::size_t b = 0; b < labels.size(); ++b)
          me->setBinLabel(static_cast<int>(b + 1), labels[b]);
        return me;
      };
      f.ladderTrackFound = bookLadder("ladder_trackFound", kLadderTrackFoundLabels);
      f.ladderTrackMissing = bookLadder("ladder_trackMissing", kLadderTrackMissingLabels);
      f.responseTrackFound = booker.book1D("response_trackFound", "candidate / truth energy;E_{cand}/E_{truth};", 60, 0., 3.);
      f.responseTrackMissing =
          booker.book1D("response_trackMissing", "candidate / truth energy;E_{cand}/E_{truth};", 60, 0., 3.);
      f.foreignFraction = booker.book1D("foreign_fraction", "foreign energy share of the candidate;fraction;", 50, 0., 1.);
    }
  }
  std::vector<std::string> types = kPFTypes;
  types.push_back("all");
  for (auto const& t : types) {
    for (auto const& region : regions) {
      booker.setCurrentFolder(dirName_ + "/reco/" + t + "/" + region);
      RecoFolder& f = reco_[t][region];
      f.reco = bookAxes(booker, "num_reco", false);
      for (auto const& c : kClasses)
        f.classes.push_back(bookAxes(booker, "num_" + c, false));
      f.ownFraction = booker.book1D("own_fraction", "energy share of the level branch;fraction;", 50, 0., 1.);
    }
  }
}

void PFCandidateTruthValidator::fillTruth(TruthFolder& f, truth::PFCandidateTruthRecord const& r) const {
  f.simul.fill(r.pt, r.eta, r.caloEta);
  for (std::size_t k = 0; k < kRungs.size(); ++k) {
    auto const& rung = kRungs[k];
    const bool applies = rung.alwaysExpected ? true : rungApplies(r, rung.passed);
    if (!applies)
      continue;
    f.expected[k].fill(r.pt, r.eta, r.caloEta);
    if (r.has(rung.passed))
      f.passed[k].fill(r.pt, r.eta, r.caloEta);
  }
  // Cut-flows. A rung that does not apply is skipped; the first failing rung ends it.
  auto ladder = [&](MonitorElement* me, std::vector<PFRung> const& rungs) {
    me->Fill(0.);
    for (std::size_t k = 0; k < rungs.size(); ++k) {
      if (!rungApplies(r, rungs[k]))
        continue;
      if (!r.has(rungs[k]))
        break;
      me->Fill(static_cast<double>(k + 1));
    }
  };
  const bool trackBranch = r.has(PFRung::TrackExpected) && r.has(PFRung::TrackFound);
  const bool missingBranch = r.has(PFRung::TrackExpected) && !r.has(PFRung::TrackFound);
  if (trackBranch)
    ladder(f.ladderTrackFound, kLadderTrackFound);
  else
    ladder(f.ladderTrackMissing, kLadderTrackMissing);  // neutral particles share this branch
  if (r.has(PFRung::CandidateFound)) {
    (trackBranch ? f.responseTrackFound : f.responseTrackMissing)->Fill(r.energyRatio);
    f.foreignFraction->Fill(r.foreignFraction);
  }
  (void)missingBranch;
}

void PFCandidateTruthValidator::fillReco(RecoFolder& f, truth::PFCandidateRecord const& c) const {
  f.reco.fill(c.pt, c.eta, c.eta);
  if (c.candidateClass < kClasses.size())
    f.classes[c.candidateClass].fill(c.pt, c.eta, c.eta);
  f.ownFraction->Fill(c.ownFraction);
}

void PFCandidateTruthValidator::analyze(edm::Event const& event, edm::EventSetup const&) {
  for (auto const& r : event.get(truthToken_)) {
    const std::string s = speciesOf(r.pdgId);
    const std::string region = r.region < kRegions.size() ? kRegions[r.region] : "forward";
    fillTruth(truth_[s][region], r);
    fillTruth(truth_[s]["all"], r);
    fillTruth(truth_["all"][region], r);
    fillTruth(truth_["all"]["all"], r);
  }
  for (auto const& c : event.get(recoToken_)) {
    const std::string t = c.pfType >= 0 && static_cast<std::size_t>(c.pfType) < kPFTypes.size() ? kPFTypes[c.pfType] : "X";
    const std::string region = c.region < kRegions.size() ? kRegions[c.region] : "forward";
    fillReco(reco_[t][region], c);
    fillReco(reco_[t]["all"], c);
    fillReco(reco_["all"][region], c);
    fillReco(reco_["all"]["all"], c);
  }
}

void PFCandidateTruthValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("truthRecords", edm::InputTag("pfCandidateTruthAssociator", "particleFlowTruthRecords"));
  desc.add<edm::InputTag>("candidateRecords", edm::InputTag("pfCandidateTruthAssociator", "particleFlowCandidateRecords"));
  desc.add<std::string>("dirName", "TruthInfo/Offline/PFCandidates/particleFlow");
  desc.add<int>("nintPt", 30);
  desc.add<double>("minPt", 1.);
  desc.add<double>("maxPt", 1000.);
  desc.add<int>("nintEta", 30);
  desc.add<double>("maxEta", 3.);
  descriptions.add("pfCandidateTruthValidator", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFCandidateTruthValidator);
