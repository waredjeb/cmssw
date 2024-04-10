import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
pfClusterHBHEOnlyAlpakaComparison = DQMEDAnalyzer("PFCaloGPUComparisonTask",
                                                    pfClusterToken_ref = cms.untracked.InputTag('particleFlowClusterHBHEOnly'),
                                                    pfClusterToken_target = cms.untracked.InputTag('legacyPFClusterProducerHBHEOnly'),
                                                    pfCaloGPUCompDir = cms.untracked.string("pfClusterHBHEAlpakaV")
)

pfClusterHBHEAlpakaComparison = DQMEDAnalyzer("PFCaloGPUComparisonTask",
                                                    pfClusterToken_ref = cms.untracked.InputTag('particleFlowClusterHBHE'),
                                                    pfClusterToken_target = cms.untracked.InputTag('legacyPFClusterProducer'),
                                                    pfCaloGPUCompDir = cms.untracked.string("pfClusterHBHEAlpakaV")
)

HLTpfClusterHBHEAlpakaComparisonGPUvsCPU = DQMEDAnalyzer("PFCaloGPUComparisonTask",
                                                    pfClusterToken_ref = cms.untracked.InputTag('hltParticleFlowClusterHCAL'),
                                                    pfClusterToken_target = cms.untracked.InputTag('hltParticleFlowClusterHCALCPUOnly'),
                                                    pfCaloGPUCompDir = cms.untracked.string("pfClusterHBHEAlpakaGPUvsCPU"),
                                                    GPU = cms.untracked.string("AlpakaGPU"),
                                                    CPU = cms.untracked.string("AlpakaCPU")
)
HLTpfClusterHBHEAlpakaComparisonCPUvsLegacy = DQMEDAnalyzer("PFCaloGPUComparisonTask",
                                                    pfClusterToken_ref = cms.untracked.InputTag('hltParticleFlowClusterHCALCPUOnly'),
                                                    pfClusterToken_target = cms.untracked.InputTag('hltParticleFlowClusterHCALLegacy'),
                                                    pfCaloGPUCompDir = cms.untracked.string("pfClusterHBHEAlpakaCPUvsLegacy"),
                                                    GPU = cms.untracked.string("AlpakaCPU"),
                                                    CPU = cms.untracked.string("Legacy")
)
HLTpfClusterHBHEAlpakaComparisonGPUvsLegacy = DQMEDAnalyzer("PFCaloGPUComparisonTask",
                                                    pfClusterToken_ref = cms.untracked.InputTag('hltParticleFlowClusterHCALCPUOnly'),
                                                    pfClusterToken_target = cms.untracked.InputTag('hltParticleFlowClusterHCALLegacy'),
                                                    pfCaloGPUCompDir = cms.untracked.string("pfClusterHBHEAlpakaGPUvsLegacy"),
                                                    GPU = cms.untracked.string("AlpakaGPU"),
                                                    CPU = cms.untracked.string("Legacy")
)
