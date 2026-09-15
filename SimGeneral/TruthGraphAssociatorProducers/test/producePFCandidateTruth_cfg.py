# Runs the PFCandidate truth association on a Phase-2 step3: rebuilds the truth-graph
# chain with a selection preset, the branch targets, the constituent associators (tracks,
# PFClusters, tracksters) and PFCandidateTruthAssociator, and writes the maps and the
# records to an EDM file.
#
#   cmsRun producePFCandidateTruth_cfg.py step3.root -n 10 -g D110 -p top -o pfcandTruth.root

import FWCore.ParameterSet.Config as cms
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("inputFile", nargs='?', default="step3.root", metavar='FILE')
parser.add_argument('-n', "--maxevts", type=int, default=10)
parser.add_argument('-g', "--geometry", default="D110")
parser.add_argument('-p', "--preset", default="top",
                    help="selection template naming the signal seeds; reconstructableFromSignal is empty without one")
parser.add_argument('-o', "--out", default="pfcandTruth.root")
args = parser.parse_args()
if ':' not in args.inputFile:
    args.inputFile = 'file:' + args.inputFile

process = cms.Process("PFCANDTRUTH")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.load("Configuration.Geometry.GeometryExtendedRun4%sReco_cff" % args.geometry)
process.trackerGeometry.applyAlignment = cms.bool(False)
process.load("Validation.Configuration.truthPrevalidation_cff")
process.load("SimGeneral.TruthGraphAssociatorProducers.truthGraphAssociators_cff")
from SimGeneral.TruthGraphAssociatorProducers.pfCandidateTruthAssociator_cfi import pfCandidateTruthAssociator
process.pfCandidateTruthAssociator = pfCandidateTruthAssociator.clone()

if args.preset:
    from PhysicsTools.TruthInfo.truthGraphSelections import postProcessingPSet, seedPdgIdsForPreset
    process.truthLogicalGraphProducer.postProcessing = postProcessingPSet(template=args.preset)
    process.truthBranchTargets.signalSeedPdgIds = cms.vint32(*seedPdgIdsForPreset(template=args.preset))

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(args.maxevts))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(args.inputFile))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False))

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName=cms.untracked.string(args.out),
    outputCommands=cms.untracked.vstring(
        "drop *",
        "keep *_pfCandidateTruthAssociator_*_*",
        "keep *_truthBranchTargets_*_*",
        "keep *_truthLogicalGraphProducer_*_*",
        "keep *_particleFlow_*_*",
    ),
)
process.p = cms.Path(
    process.truthGraphProducer
    + process.truthLogicalGraphProducer
    + process.detIdToRecHitMapProducer
    + process.truthLogicalGraphHitIndexProducer
    + process.truthBranchTargets
    + process.allTrackToTruthBranchAssociators
    + process.truthBranchPFClusterEcalAssociators
    + process.truthBranchPFClusterHcalAssociators
    + process.truthBranchTracksterAssociators
    + process.pfCandidateTruthAssociator
)
process.e = cms.EndPath(process.out)
