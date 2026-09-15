# PFCandidate truth validation on a Phase-2 step3: the truth-graph chain with a
# selection preset, the branch targets, the constituent associators, the PFCandidate
# associator and its DQM analyzer. Output: a DQMIO file for harvestPFCandidateTruth_cfg.py.
#
#   cmsRun validatePFCandidateTruth_cfg.py step3.root -n 100 -g D110 -p top -o pfcand_dqm.root   (gun: -p gun -l reconstructableFinalState)

import FWCore.ParameterSet.Config as cms
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("inputFile", nargs='?', default="step3.root", metavar='FILE')
parser.add_argument('-n', "--maxevts", type=int, default=10)
parser.add_argument('-g', "--geometry", default="D110")
parser.add_argument('-p', "--preset", default="top", help="selection template naming the signal seeds (top, gun, ...)")
parser.add_argument('-l', "--level", default="reconstructableFromSignal",
                    help="truth level of the ladder denominator; a gun has no resonance, use reconstructableFinalState")
parser.add_argument('-o', "--out", default="pfcand_dqm.root")
args = parser.parse_args()
if ':' not in args.inputFile:
    args.inputFile = 'file:' + args.inputFile

process = cms.Process("PFCANDDQM")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.load("Configuration.Geometry.GeometryExtendedRun4%sReco_cff" % args.geometry)
process.load("DQMServices.Core.DQMStore_cfi")
process.trackerGeometry.applyAlignment = cms.bool(False)
process.load("Validation.Configuration.truthPrevalidation_cff")
process.load("SimGeneral.TruthGraphAssociatorProducers.truthGraphAssociators_cff")
process.load("Validation.TruthInfo.pfCandidateTruthValidation_cff")

if args.preset:
    from PhysicsTools.TruthInfo.truthGraphSelections import postProcessingPSet, seedPdgIdsForPreset
    process.truthLogicalGraphProducer.postProcessing = postProcessingPSet(template=args.preset)
    process.truthBranchTargets.signalSeedPdgIds = cms.vint32(*seedPdgIdsForPreset(template=args.preset))

process.pfCandidateTruthAssociator.targets = cms.InputTag(
    "truthBranchTargets", "truthToRecoTargets" + args.level[0].upper() + args.level[1:])

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(args.maxevts))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(args.inputFile))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False))
process.dqmOut = cms.OutputModule("DQMRootOutputModule", fileName=cms.untracked.string(args.out))

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
    + process.pfCandidateTruthValidationSequence
)
process.e = cms.EndPath(process.dqmOut)
