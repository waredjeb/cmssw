import FWCore.ParameterSet.Config as cms

process = cms.Process('ParticleFlowDQMOffline')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# load DQM
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")

# load jet correctors
process.load('JetMETCorrections.Configuration.JetCorrectors_cff')

# my analyzer
process.load('DQMOffline.ParticleFlow.runBasic_cfi')

# Setup Global Tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        connectionTimeout = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    DumpStat = cms.untracked.bool(False),
    JsonDumpFileName = cms.untracked.string(''),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    appendToDataLabel = cms.string(''),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    frontierKey = cms.untracked.string(''),
    globaltag = cms.string('150X_mcRun4_realistic_v1'),
    pfnPostfix = cms.untracked.string(''),
    pfnPrefix = cms.untracked.string(''),
    recordsToDebug = cms.untracked.vstring(),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet(
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            record = cms.string('SiPixelGenErrorDBObjectRcd'),
            snapshotTime = cms.string('2023-05-16 20:00:00'),
            tag = cms.string('SiPixelGenErrorDBObject_phase2_IT_v7.1.1_25x100_v1_mc')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            record = cms.string('SiPixelLorentzAngleRcd'),
            snapshotTime = cms.string('2023-05-16 20:00:00.000'),
            tag = cms.string('SiPixelLorentzAngle_phase2_IT_v7.1.1_25x100_v1_mc')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('forWidth'),
            record = cms.string('SiPixelLorentzAngleRcd'),
            snapshotTime = cms.string('2023-12-02 15:55:00.000'),
            tag = cms.string('SiPixelLorentzAngle_phase2_IT_v7.1.1_25x100_empty_mc')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('fromAlignment'),
            record = cms.string('SiPixelLorentzAngleRcd'),
            snapshotTime = cms.string('2023-12-02 15:55:00.000'),
            tag = cms.string('SiPixelLorentzAngle_phase2_IT_v7.1.1_25x100_empty_mc')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            record = cms.string('SiPixelLorentzAngleSimRcd'),
            snapshotTime = cms.string('2024-03-07 21:00:00.000'),
            tag = cms.string('SiPixelSimLorentzAngle_phase2_IT_v7.1.1_25x100_v1_mc')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            record = cms.string('SiPixelTemplateDBObjectRcd'),
            snapshotTime = cms.string('2023-05-16 20:00:00'),
            tag = cms.string('SiPixelTemplateDBObject_phase2_IT_v7.1.1_25x100_v1_mc')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            record = cms.string('TrackerAlignmentErrorExtendedRcd'),
            snapshotTime = cms.string('2024-09-12 15:37:00'),
            tag = cms.string('TrackerAlignmentErrorsExtended_Upgrade2026_T33_design_v1')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            record = cms.string('TrackerAlignmentRcd'),
            snapshotTime = cms.string('2024-09-12 15:37:00'),
            tag = cms.string('TrackerAlignment_Upgrade2026_T33_design_v1')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            record = cms.string('TrackerSurfaceDeformationRcd'),
            snapshotTime = cms.string('2023-03-16 15:30:00'),
            tag = cms.string('TrackerSurfaceDeformations_Upgrade2026_Zero')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            record = cms.string('EcalIntercalibConstantsRcd'),
            tag = cms.string('EcalIntercalibConstants_TL1000_upgrade_8deg_v2_mc')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            record = cms.string('EcalIntercalibConstantsMCRcd'),
            tag = cms.string('EcalIntercalibConstantsMC_TL1000_upgrade_8deg_v2_mc')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            record = cms.string('EcalLaserAPDPNRatiosRcd'),
            tag = cms.string('EcalLaserAPDPNRatios_TL1000_upgrade_8deg_mc')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            record = cms.string('EcalPedestalsRcd'),
            tag = cms.string('EcalPedestals_TL1000_upgradeTIA_8deg_mc')
        ),
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            record = cms.string('EcalTPGLinearizationConstRcd'),
            tag = cms.string('EcalTPGLinearizationConst_TL1000_upgrade_8deg_mc')
        ), 
        template = cms.PSetTemplate(
            connect = cms.string(''),
            label = cms.untracked.string(''),
            record = cms.string(''),
            refreshTime = cms.uint64(18446744073709551615),
            snapshotTime = cms.string(''),
            tag = cms.string('')
        )
    )
)

# Here we explicitly override the jet energy corrections (JECs) in a Global Tag
process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(
    record = cms.string("JetCorrectionsRecord"),
    tag = cms.string("JetCorrectorParametersCollection_Winter25Prompt25_RunC_V1_DATA_AK4PFPuppi_v1"),
    label = cms.untracked.string('AK4PFPuppi'),
    connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
  )
)

with open('fileList.log') as f:
    lines = f.readlines()
#Input source
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(lines))

# "CorrectedPFJetProducer" module applies the jet energy
# corrections on the jet collection and sort the collection
# according to pt
process.hltFixedGridRhoFastjetAll = cms.EDProducer("FixedGridRhoProducerFastjet",
    gridSpacing = cms.double(0.55),
    maxRapidity = cms.double(5.0),
    pfCandidatesTag = cms.InputTag("hltParticleFlowTmp")
)
process.hltAK4PFJetCorrectorL1 = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("hltFixedGridRhoFastjetAll")
)


process.hltAK4PFJetCorrectorL2 = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2Relative')
)


process.hltAK4PFJetCorrectorL3 = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L3Absolute')
)
process.hltAK4PFJetCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("hltAK4PFJetCorrectorL1", "hltAK4PFJetCorrectorL2", "hltAK4PFJetCorrectorL3")
)
process.hltAK4PFJetsCorrected = cms.EDProducer("CorrectedPFJetProducer",
    correctors = cms.VInputTag("hltAK4PFJetCorrector"),
    src = cms.InputTag("hltAK4PFJets")
)

###################################################################
# Data certification GoldenJSON filtering
###################################################################
goldenJSONPath=""
if goldenJSONPath != "":
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = goldenJSONPath).getVLuminosityBlockRange()

from DQMOffline.ParticleFlow.runBasic_cfi import *

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
                                     fileName = cms.untracked.string("OUT_step1.root"))


process.p = cms.Path(
    #process.hltFixedGridRhoFastjetAll+
    #process.hltAK4PFJetCorrectorL1+
    #process.hltAK4PFJetCorrectorL2+
    #process.hltAK4PFJetCorrectorL3+
    #process.hltAK4PFJetCorrector+
    #process.hltAK4PFJetsCorrected+
    process.PFAnalyzer)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)


## Schedule definition
process.schedule = cms.Schedule(
    process.p,
    process.DQMoutput_step
    )









