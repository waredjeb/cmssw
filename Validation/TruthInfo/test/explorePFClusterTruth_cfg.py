# Driver for TruthGraphPFClusterExplorer on a Phase-2 step3: rebuilds the truth-graph
# chain, the branch targets and the PFCluster associators, then runs the explorer.
#
#   cmsRun explorePFClusterTruth_cfg.py path/to/step3.root -n 5 -g D110

import FWCore.ParameterSet.Config as cms
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("inputFile", nargs='?', default="step3.root", metavar='FILE')
parser.add_argument('-n', "--maxevts", type=int, default=5)
parser.add_argument('-g', "--geometry", default="D110", help="Run4 geometry tag of the sample")
args = parser.parse_args()
if '/' not in args.inputFile and ':' not in args.inputFile:
    args.inputFile = 'file:' + args.inputFile
elif ':' not in args.inputFile:
    args.inputFile = 'file:' + args.inputFile

process = cms.Process("PFCLEXPLORE")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.threshold = "INFO"
process.MessageLogger.cerr.TruthGraphPFClusterExplorer = cms.untracked.PSet(limit=cms.untracked.int32(-1))
process.load("Configuration.Geometry.GeometryExtendedRun4%sReco_cff" % args.geometry)
process.trackerGeometry.applyAlignment = cms.bool(False)

# truthGraphProducer -> truthLogicalGraphProducer -> detIdToRecHitMapProducer -> hit index
process.load("Validation.Configuration.truthPrevalidation_cff")
# truthBranchTargets + the associators (tracks, vertices, tracksters, PFClusters)
process.load("SimGeneral.TruthGraphAssociatorProducers.truthGraphAssociators_cff")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(args.maxevts))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(args.inputFile))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False))

from Validation.TruthInfo.truthGraphPFClusterExplorer_cfi import truthGraphPFClusterExplorer
process.explorer = truthGraphPFClusterExplorer.clone()

process.p = cms.Path(
    process.truthGraphProducer
    + process.truthLogicalGraphProducer
    + process.detIdToRecHitMapProducer
    + process.truthLogicalGraphHitIndexProducer
    + process.truthBranchTargets
    + process.truthBranchPFClusterEcalAssociators
    + process.truthBranchPFClusterHcalAssociators
    + process.explorer
)
