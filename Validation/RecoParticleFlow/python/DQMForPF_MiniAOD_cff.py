import FWCore.ParameterSet.Config as cms

from Validation.RecoParticleFlow.particleFlowDQM_cff import pfJetAnalyzerDQM,pfJetAnalyzerHLTDQM 
#from Validation.RecoParticleFlow.particleFlowDQM_cff import pfJetAnalyzerDQM
from Validation.RecoParticleFlow.particleFlowDQM_cff import pfPuppiJetAnalyzerDQM
from Validation.RecoParticleFlow.particleFlowDQM_cff import pfJetDQMPostProcessor
from Validation.RecoParticleFlow.particleFlowDQM_cff import PFCandAnalyzerDQM
from Validation.RecoParticleFlow.offsetAnalyzerDQM_cff import offsetAnalyzerDQM
from Validation.RecoParticleFlow.particleFlowDQM_cff import PFCandAnalyzerHLTDQM
from Validation.RecoParticleFlow.offsetAnalyzerDQM_cff import offsetDQMPostProcessor
from Validation.RecoParticleFlow.pfCaloGPUComparisonTask_cfi import HLTpfClusterHBHEAlpakaComparisonGPUvsCPU, HLTpfClusterHBHEAlpakaComparisonCPUvsLegacy, HLTpfClusterHBHEAlpakaComparisonGPUvsLegacy
from Validation.RecoParticleFlow.particleFlowDQM_cff import pfJetHLTDQMPostProcessor
from RecoLocalCalo.HcalRecAlgos.hcalChannelPropertiesESProd_cfi import hcalChannelPropertiesESProd
# Use also other POGs' analyzers for extended checks
from Validation.RecoMET.METRelValForDQM_cff import *
from Validation.RecoJets.JetValidation_cff import *

pfBuilder = cms.PSet(
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

hltParticleFlowRecHitHBHE = cms.EDProducer("LegacyPFRecHitProducer",
        src = cms.InputTag("hltParticleFlowRecHitHBHESoACPUSerial")
  )

hltParticleFlowClusterHBHE = cms.EDProducer("LegacyPFClusterProducer",
    src = cms.InputTag("hltParticleFlowClusterHBHESoACPUSerial"),
    pfClusterBuilder = pfBuilder, 
    usePFThresholdsFromDB = cms.bool(True),
    recHitsSource = cms.InputTag("hltParticleFlowRecHitHBHE"),
    PFRecHitsLabelIn = cms.InputTag("hltParticleFlowRecHitHBHESoACPUSerial")
)


hltParticleFlowClusterHCALCPUOnly = cms.EDProducer("PFMultiDepthClusterProducer",
    clustersSource = cms.InputTag("hltParticleFlowClusterHBHE"),
    energyCorrector = cms.PSet(

    ),
    pfClusterBuilder = cms.PSet(
        algoName = cms.string('PFMultiDepthClusterizer'),
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
        minFractionToKeep = cms.double(1e-07),
        nSigmaEta = cms.double(2.0),
        nSigmaPhi = cms.double(2.0)
    ),
    positionReCalc = cms.PSet(

    ),
    usePFThresholdsFromDB = cms.bool(True)
)

HLTpfClusterHBHEAlpakaComparisonGPUvsCPU = DQMEDAnalyzer("PFCaloGPUComparisonTask",
                                                    pfClusterToken_ref = cms.untracked.InputTag('hltParticleFlowClusterHCAL'),
                                                    pfClusterToken_target = cms.untracked.InputTag('hltParticleFlowClusterHCALCPUOnly'),
                                                    pfCaloGPUCompDir = cms.untracked.string("pfClusterHBHEAlpakaGPUvsCPU"),
                                                    GPU = cms.untracked.string("AlpakaGPU"),
                                                    CPU = cms.untracked.string("AlpakaCPU")
)

DQMHLTPF = cms.Sequence(
    hltParticleFlowRecHitHBHE+
    hltParticleFlowClusterHBHE+
    hltParticleFlowClusterHCALCPUOnly+
    HLTpfClusterHBHEAlpakaComparisonGPUvsCPU
)


#DQMHLTPF = cms.Sequence(
#    #pfJetAnalyzerHLTDQM+
#    #PFCandAnalyzerHLTDQM+
#    hltParticleFlowRecHitHBHE+
#    hltParticleFlowClusterHBHE 
##    hltParticleFlowClusterHBHE+
##    hltParticleFlowClusterHBHE+
#    #HLTpfClusterHBHEAlpakaComparisonCPUvsLegacy+
#    #HLTpfClusterHBHEAlpakaComparisonGPUvsLegacy
#)

DQMHarvestHLTPF = cms.Sequence(
  pfJetHLTDQMPostProcessor
  )

