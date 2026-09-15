# DQM analyzer and harvester of the PFCandidate truth validation. The analyzer counts
# from the records PFCandidateTruthAssociator persists; the harvester turns every rung
# into an efficiency against its own denominator and every candidate class into a
# fraction of the candidates.

import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester
from Validation.TruthInfo.pfCandidateTruthValidator_cfi import pfCandidateTruthValidator

_dir = pfCandidateTruthValidator.dirName.value()
_regions = ["barrel", "transition", "endcap", "forward", "all"]
_species = ["pi", "K", "p", "e", "mu", "gamma", "pi0", "n", "K0L", "other", "all"]
_pftypes = ["X", "h", "e", "mu", "gamma", "h0", "h_HF", "egamma_HF", "all"]
_classes = ["unmatched", "matched", "merged", "split", "other", "unevaluated"]
_rungs = ["trackFound", "trackInCandidate", "ecalCollected", "ecalLinked", "hcalCollected", "hcalLinked",
          "hgcalCollected", "hgcalLinked", "caloSameCandidate", "candidateFound", "merged", "partial",
          "clean", "pdgCorrect"]
_axes = ["pt", "eta", "caloeta"]


def _efficiencyStrings():
    out = []
    for rung in _rungs:
        for axis in _axes:
            out.append("eff_%s_vs_%s 'efficiency of %s vs %s' num_%s_%s num_%s_expected_%s"
                       % (rung, axis, rung, axis, rung, axis, rung, axis))
    return out


def _classStrings():
    out = []
    for cls in _classes:
        for axis in ("pt", "eta"):
            out.append("frac_%s_vs_%s 'fraction of %s candidates vs %s' num_%s_%s num_reco_%s"
                       % (cls, axis, cls, axis, cls, axis, axis))
    return out


pfCandidateTruthValidationSequence = cms.Sequence(pfCandidateTruthValidator)

pfCandidateTruthPostProcessor = DQMEDHarvester(
    "DQMGenericClient",
    subDirs=cms.untracked.vstring(*["%s/truth/%s/%s" % (_dir, s, r) for s in _species for r in _regions]),
    efficiency=cms.vstring(*_efficiencyStrings()),
    resolution=cms.vstring(),
    verbose=cms.untracked.uint32(0),
    outputFileName=cms.untracked.string(""),
)
pfCandidateRecoPostProcessor = DQMEDHarvester(
    "DQMGenericClient",
    subDirs=cms.untracked.vstring(*["%s/reco/%s/%s" % (_dir, t, r) for t in _pftypes for r in _regions]),
    efficiency=cms.vstring(*_classStrings()),
    resolution=cms.vstring(),
    verbose=cms.untracked.uint32(0),
    outputFileName=cms.untracked.string(""),
)
pfCandidateTruthHarvestingSequence = cms.Sequence(pfCandidateTruthPostProcessor + pfCandidateRecoPostProcessor)
