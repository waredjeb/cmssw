import FWCore.ParameterSet.Config as cms

from Validation.RecoParticleFlow.particleFlowDQM_cff import pfJetAnalyzerDQM,pfJetAnalyzerHLTDQM 
#from Validation.RecoParticleFlow.particleFlowDQM_cff import pfJetAnalyzerDQM
from Validation.RecoParticleFlow.particleFlowDQM_cff import pfPuppiJetAnalyzerDQM
from Validation.RecoParticleFlow.particleFlowDQM_cff import pfJetDQMPostProcessor
from Validation.RecoParticleFlow.particleFlowDQM_cff import PFCandAnalyzerDQM
from Validation.RecoParticleFlow.offsetAnalyzerDQM_cff import offsetAnalyzerDQM
from Validation.RecoParticleFlow.particleFlowDQM_cff import PFCandAnalyzerHLTDQM
from Validation.RecoParticleFlow.offsetAnalyzerDQM_cff import offsetDQMPostProcessor
#from DQM.ParticleFlow.pfCaloGPUComparisonTask_cfi import HLTpfClusterHBHEAlpakaComparisonGPUvsCPU#, HLTpfClusterHBHEAlpakaComparisonCPUvsLegacy, HLTpfClusterHBHEAlpakaComparisonGPUvsLegacy
from Validation.RecoParticleFlow.particleFlowDQM_cff import pfJetHLTDQMPostProcessor
# Use also other POGs' analyzers for extended checks
from Validation.RecoMET.METRelValForDQM_cff import *
from Validation.RecoJets.JetValidation_cff import *


DQMOfflinePF = cms.Sequence(
  pfJetAnalyzerDQM +
  pfPuppiJetAnalyzerDQM +
  offsetAnalyzerDQM +
  PFCandAnalyzerDQM
)

DQMHarvestPF = cms.Sequence(
  pfJetDQMPostProcessor +
  offsetDQMPostProcessor
)

# MET & Jets sequence
DQMOfflinePFExtended = cms.Sequence(
    METValidationMiniAOD +
    JetValidationMiniAOD
)


DQMHLTPF = cms.Sequence(
#    pfJetAnalyzerHLTDQM+
    PFCandAnalyzerHLTDQM
    #HLTpfClusterHBHEAlpakaComparisonGPUvsCPU
    
)

DQMHarvestHLTPF = cms.Sequence(
  pfJetHLTDQMPostProcessor
  )

