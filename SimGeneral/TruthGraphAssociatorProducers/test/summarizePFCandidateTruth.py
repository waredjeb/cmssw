#!/usr/bin/env python3
"""FWLite summary of the PFCandidate truth records: ladder pass rates per species and
region, and the candidate classes per PF type."""
import sys
from collections import Counter, defaultdict
import ROOT
from DataFormats.FWLite import Events, Handle

RUNGS = ["TrackExpected", "TrackFound", "TrackInCandidate", "EcalExpected", "EcalCollected", "EcalLinked",
         "HcalExpected", "HcalCollected", "HcalLinked", "HgcalExpected", "HgcalCollected", "HgcalLinked",
         "CaloSameCandidate", "CandidateFound", "Merged", "Clean", "PdgEvaluated", "PdgCorrect"]
BIT = {n: 1 << i for i, n in enumerate(RUNGS)}
REGION = ["barrel", "transition", "endcap", "forward"]
CLASS = ["Unmatched", "Matched", "Merged", "Split", "Other", "Unevaluated"]
PFTYPE = ["X", "h", "e", "mu", "gamma", "h0", "h_HF", "egamma_HF"]
SPECIES = {211: "pi", 321: "K", 2212: "p", 11: "e", 13: "mu", 22: "gamma", 111: "pi0", 2112: "n", 130: "K0L"}

events = Events(sys.argv[1])
truthH, candH = Handle("std::vector<truth::PFCandidateTruthRecord>"), Handle("std::vector<truth::PFCandidateRecord>")
label = ("pfCandidateTruthAssociator", "particleFlowTruthRecords")
clabel = ("pfCandidateTruthAssociator", "particleFlowCandidateRecords")
ladder = defaultdict(Counter)   # (species, region) -> rung counts
classes = defaultdict(Counter)  # pftype -> class counts
n = 0
for ev in events:
    n += 1
    ev.getByLabel(label, truthH); ev.getByLabel(clabel, candH)
    for r in truthH.product():
        key = (SPECIES.get(abs(r.pdgId), "other"), REGION[r.region])
        ladder[key]["all"] += 1
        for name, bit in BIT.items():
            if r.rungs & bit: ladder[key][name] += 1
    for c in candH.product():
        classes[PFTYPE[c.pfType]][CLASS[c.candidateClass]] += 1
print("events:", n)
print("\n%-8s %-10s %6s | %5s %5s %5s | %5s %5s %5s | %5s %5s %5s | %5s %5s %5s | %5s %5s %5s | %5s %5s %5s" % (
    "species", "region", "n", "trkE", "trkF", "trkC", "ecE", "ecC", "ecL", "hcE", "hcC", "hcL", "hgE", "hgC", "hgL", "same", "cand", "merg", "clean", "pdgE", "pdgOK"))
for key in sorted(ladder, key=lambda k: (-ladder[k]["all"], k)):
    c = ladder[key]
    if c["all"] < 20: continue
    print("%-8s %-10s %6d | %5d %5d %5d | %5d %5d %5d | %5d %5d %5d | %5d %5d %5d | %5d %5d %5d | %5d %5d %5d" % (
        key[0], key[1], c["all"], c["TrackExpected"], c["TrackFound"], c["TrackInCandidate"],
        c["EcalExpected"], c["EcalCollected"], c["EcalLinked"], c["HcalExpected"], c["HcalCollected"], c["HcalLinked"],
        c["HgcalExpected"], c["HgcalCollected"], c["HgcalLinked"], c["CaloSameCandidate"], c["CandidateFound"], c["Merged"],
        c["Clean"], c["PdgEvaluated"], c["PdgCorrect"]))
print("\ncandidate classes per PF type:")
for t in PFTYPE:
    if classes[t]: print("  %-6s %s" % (t, dict(classes[t])))
