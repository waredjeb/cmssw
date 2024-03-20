import FWCore.ParameterSet.Config as cms

from Validation.RecoParticleFlow.particleFlowDQM_cff import pfJetAnalyzerDQM,pfJetAnalyzerHLTDQM 
#from Validation.RecoParticleFlow.particleFlowDQM_cff import pfJetAnalyzerDQM
from Validation.RecoParticleFlow.particleFlowDQM_cff import pfPuppiJetAnalyzerDQM
from Validation.RecoParticleFlow.particleFlowDQM_cff import pfJetDQMPostProcessor
from Validation.RecoParticleFlow.particleFlowDQM_cff import PFCandAnalyzerDQM
from Validation.RecoParticleFlow.offsetAnalyzerDQM_cff import offsetAnalyzerDQM
from Validation.RecoParticleFlow.particleFlowDQM_cff import PFCandAnalyzerHLTDQM
from Validation.RecoParticleFlow.offsetAnalyzerDQM_cff import offsetDQMPostProcessor
from DQM.ParticleFlow.pfCaloGPUComparisonTask_cfi import HLTpfClusterHBHEAlpakaComparisonGPUvsCPU, HLTpfClusterHBHEAlpakaComparisonCPUvsLegacy, HLTpfClusterHBHEAlpakaComparisonGPUvsLegacy
from Validation.RecoParticleFlow.particleFlowDQM_cff import pfJetHLTDQMPostProcessor
# Use also other POGs' analyzers for extended checks
from Validation.RecoMET.METRelValForDQM_cff import *
from Validation.RecoJets.JetValidation_cff import *

pfClusterBuilder = cms.PSet(                                                                                                                                                                                                                                                                                                                  
    algoName = cms.string('Basic2DGenericPFlowClusterizer'),                                                                                                                                                                                                                                                               
    allCellsPositionCalc = cms.PSet(                                                                                                                                                                                                                                                                                       
        algoName = cms.string('Basic2DGenericPFlowPositionCalc'),                                                                                                                                                                                                                                                          
        logWeightDenominatorByDetector = cms.VPSet(                                                                                                                                                                                                                                                                        
            cms.PSet(                                                                                                                                                                                                                                                                                                      
                depths = cms.vint32(1, 2, 3, 4),                                                                                                                                                                                                                                                                           
                detector = cms.string('HCAL_BARREL1'),                                                                                                                                                                                                                                                                     
                logWeightDenominator = cms.vdouble(0.4, 0.3, 0.3, 0.3)                                                                                                                                                                                                                                                     
            ),                                                                                                                                                                                                                                                                                                             
            cms.PSet(                                                                                                                                                                                                                                                                                                      
                depths = cms.vint32(                                                                                                                                                                                                                                                                                       
                    1, 2, 3, 4, 5,                                                                                                                                                                                                                                                                                         
                    6, 7                                                                                                                                                                                                                                                                                                   
                ),                                                                                                                                                                                                                                                                                                         
                detector = cms.string('HCAL_ENDCAP'),                                                                                                                                                                                                                                                                      
                logWeightDenominator = cms.vdouble(                                                                                                                                                                                                                                                                        
                    0.1, 0.2, 0.2, 0.2, 0.2,                                                                                                                                                                                                                                                                               
                    0.2, 0.2                                                                                                                                                                                                                                                                                               
                )                                                                                                                                                                                                                                                                                                          
            )                                                                                                                                                                                                                                                                                                              
        ),
        minAllowedNormalization = cms.double(1e-09),
        minFractionInCalc = cms.double(1e-09),
        posCalcNCrystals = cms.int32(-1)
    ),
    clusterTimeResFromSeed = cms.bool(False),
    excludeOtherSeeds = cms.bool(True),
    maxIterations = cms.uint32(5),
    maxNSigmaTime = cms.double(10.0),
    minChi2Prob = cms.double(0.0),
    minFracTot = cms.double(1e-20),
    minFractionToKeep = cms.double(1e-07),
    positionCalc = cms.PSet(
        algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
        logWeightDenominatorByDetector = cms.VPSet(
            cms.PSet(
                depths = cms.vint32(1, 2, 3, 4),
                detector = cms.string('HCAL_BARREL1'),
                logWeightDenominator = cms.vdouble(0.4, 0.3, 0.3, 0.3)
            ),
            cms.PSet(
                depths = cms.vint32(
                    1, 2, 3, 4, 5,
                    6, 7
                ),
                detector = cms.string('HCAL_ENDCAP'),
                logWeightDenominator = cms.vdouble(
                    0.1, 0.2, 0.2, 0.2, 0.2,
                    0.2, 0.2
                )
            )
        ),
        minAllowedNormalization = cms.double(1e-09),
        minFractionInCalc = cms.double(1e-09),
        posCalcNCrystals = cms.int32(5)
    ),
    recHitEnergyNorms = cms.VPSet(
        cms.PSet(
            depths = cms.vint32(1, 2, 3, 4),
            detector = cms.string('HCAL_BARREL1'),
            recHitEnergyNorm = cms.vdouble(0.4, 0.3, 0.3, 0.3)
        ),
        cms.PSet(
            depths = cms.vint32(
                1, 2, 3, 4, 5,
                6, 7
            ),
            detector = cms.string('HCAL_ENDCAP'),
            recHitEnergyNorm = cms.vdouble(
                0.1, 0.2, 0.2, 0.2, 0.2,
                0.2, 0.2
            )
        )
    ),
    showerSigma = cms.double(10.0),
    stoppingTolerance = cms.double(1e-08),
    timeResolutionCalcBarrel = cms.PSet(
        constantTerm = cms.double(2.82),
        constantTermLowE = cms.double(4.24),
        corrTermLowE = cms.double(0.0),
        noiseTerm = cms.double(21.86),
        noiseTermLowE = cms.double(8.0),
        threshHighE = cms.double(15.0),
        threshLowE = cms.double(6.0)
    ),
    timeResolutionCalcEndcap = cms.PSet(
        constantTerm = cms.double(2.82),
        constantTermLowE = cms.double(4.24),
        corrTermLowE = cms.double(0.0),
        noiseTerm = cms.double(21.86),
        noiseTermLowE = cms.double(8.0),
        threshHighE = cms.double(15.0),
        threshLowE = cms.double(6.0)
    ),
    timeSigmaEB = cms.double(10.0),
    timeSigmaEE = cms.double(10.0)
)



hltParticleFlowClusterHBHE = cms.EDProducer("LegacyPFClusterProducer",
    src = cms.InputTag("hltParticleFlowClusterHBHESoA"),
    pfClusterBuilder = pfClusterBuilder, 
    usePFThresholdsFromDB = cms.bool(True),
    recHitsSource = cms.InputTag("hltParticleFlowRecHitHBHE"),
    PFRecHitsLabelIn = cms.InputTag("hltParticleFlowRecHitHBHESoA")
)

hltParticleFlowClusterHBHECPUOnly = cms.EDProducer("LegacyPFClusterProducer",
    src = cms.InputTag("hltParticleFlowClusterHBHESoACPUSerial"),
    pfClusterBuilder = pfClusterBuilder, 
    usePFThresholdsFromDB = cms.bool(True),
    recHitsSource = cms.InputTag("hltParticleFlowRecHitHBHE"),
    PFRecHitsLabelIn = cms.InputTag("hltParticleFlowRecHitHBHESoACPUSerial")
)

HLTpfClusterHBHEAlpakaComparisonGPUvsCPU= DQMEDAnalyzer("PFCaloGPUComparison",
                                                    pfClusterToken_ref = cms.untracked.InputTag('hltParticleFlowClusterHBHESoA'),
                                                    pfClusterToken_target = cms.untracked.InputTag('hltParticleFlowClusterHBHESoACPUSerial'),
                                                    pfRecHitsToken_ref = cms.untracked.InputTag('hltParticleFlowRecHitHBHESoA'),
                                                    pfRecHitsToken_target = cms.untracked.InputTag('hltParticleFlowRecHitHBHESoACPUSerial'),
                                                    pfCaloGPUCompDir = cms.untracked.string("pfClusterHBHEAlpakaV"),
                                                    GPU = cms.untracked.string("GPUAlpaka"),
                                                    CPU = cms.untracked.string("CPUAlpaka")
                                                    )

#HLTpfClusterHBHEAlpakaComparisonGPUvsCPU= DQMEDAnalyzer("PFCaloGPUComparisonTask",
#                                                    pfClusterToken_ref = cms.untracked.InputTag('hltParticleFlowClusterHCAL'),
#                                                    pfClusterToken_target = cms.untracked.InputTag('hltParticleFlowClusterHCALCPUOnly'),
#                                                    pfCaloGPUCompDir = cms.untracked.string("pfClusterHBHEAlpakaV"),
#                                                    GPU = cms.untracked.string("GPUAlpaka"),
#                                                    CPU = cms.untracked.string("CPUAlpaka")
#)

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
#    PFCandAnalyzerHLTDQM+
#    hltParticleFlowClusterHBHE+
#    hltParticleFlowClusterHBHE+
#    hltParticleFlowClusterHBHECPUOnly+
    HLTpfClusterHBHEAlpakaComparisonGPUvsCPU
    #HLTpfClusterHBHEAlpakaComparisonCPUvsLegacy+
    #HLTpfClusterHBHEAlpakaComparisonGPUvsLegacy
)

DQMHarvestHLTPF = cms.Sequence(
  pfJetHLTDQMPostProcessor
  )

