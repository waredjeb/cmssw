# Harvests the DQMIO output of validatePFCandidateTruth_cfg.py into the rung
# efficiencies and candidate-class fractions.
#
#   cmsRun harvestPFCandidateTruth_cfg.py pfcand_dqm.root

import FWCore.ParameterSet.Config as cms
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("inputFile", nargs='?', default="pfcand_dqm.root", metavar='FILE')
args = parser.parse_args()
if ':' not in args.inputFile:
    args.inputFile = 'file:' + args.inputFile

process = cms.Process("PFCANDHARVEST")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("DQMServices.Core.DQMStore_cfi")
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.load("Validation.TruthInfo.pfCandidateTruthValidation_cff")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))
process.source = cms.Source("DQMRootSource", fileNames=cms.untracked.vstring(args.inputFile))
process.options = cms.untracked.PSet(numberOfThreads=cms.untracked.uint32(1))
process.dqmSaver.convention = "Offline"
process.dqmSaver.workflow = "/TruthInfo/PFCandidates/HARVEST"
process.dqmSaver.saveByRun = cms.untracked.int32(1)
process.p = cms.Path(process.pfCandidateTruthHarvestingSequence)
process.e = cms.EndPath(process.dqmSaver)
