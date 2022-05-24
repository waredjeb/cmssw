import FWCore.ParameterSet.Config as cms
process = cms.Process("validation")

"""
start customization
"""

#Enter here the Global tags
#tag =  '94X_mc2017_realistic_v10'
#Flavour plots for MC: "all" = plots for all jets ; "dusg" = plots for d, u, s, dus, g independently ; not mandatory and any combinations are possible 
#b, c, light (dusg), non-identified (NI), PU jets plots are always produced
flavPlots = "allbcldusg"
ptRanges = cms.vdouble(30.0, 150.0, 3000.0)
#ptRanges = cms.vdouble(30.0, 70.0, 3000.0)
etaRanges = cms.vdouble(0.0, 1.4, 2.4)
#Check if jets originate from PU? option recommended (only for MC)
PUid = True

#List of taggers and taginfo to be considered (see example in: DQMOffline/RecoB/python/bTagCommon_cff.py)
from DQMOffline.RecoB.bTagCommon_cff import *
tagConfig = cms.VPSet(
        cms.PSet(
            bTagGenericAnalysisBlock,
            label = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
            folder = cms.string("CSVIVFv2-StandardRho25")    # standard CSVIVFv2 sequence
        ),
        cms.PSet(
            bTagGenericAnalysisBlock,
            label = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTagsRho90"),
            folder = cms.string("CSVIVFv2-StandardRho90")    # standard CSVIVFv2 sequence
        ),
        cms.PSet(
            bTagGenericAnalysisBlock,
            label = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTagsRho9999"),
            folder = cms.string("CSVIVFv2-StandardRho9999")    # standard CSVIVFv2 sequence
            ),
        cms.PSet(
            bTagGenericAnalysisBlock,
            label = cms.InputTag("pfDeepCSVJetTags:probb"),
            folder = cms.string("deepCSV_probb-StandardRho25")
            ),
        cms.PSet(
            bTagGenericAnalysisBlock,
            label = cms.InputTag("pfDeepCSVJetTagsRho90:probb"),
            folder = cms.string("deepCSV_probb-StandardRho90")
            ),
        cms.PSet(
            bTagGenericAnalysisBlock,
            label = cms.InputTag("pfDeepCSVJetTagsRho9999:probb"),
            folder = cms.string("deepCSV_probb-StandardRho9999")
            )
,
#         cms.PSet(
#             bTagGenericAnalysisBlock,
#             label = cms.InputTag("pfDeepFlavourJetTags:probb"),
#             folder = cms.string("deepJet_probb-StandardRho25")
#             )
# ,
#         cms.PSet(
#             bTagGenericAnalysisBlock,
#             label = cms.InputTag("pfDeepFlavourJetTagsRho90:probb"),
#             folder = cms.string("deepJet_probb-StandardRho90")
#             ),
#         cms.PSet(
#             bTagGenericAnalysisBlock,
#             label = cms.InputTag("pfDeepFlavourJetTagsRho9999:probb"),
#             folder = cms.string("deepJet_probb-StandardRho9999")
#             )
)

for i in range(0, 1) :
    tagConfig.append(
        cms.PSet(
            bTagGenericAnalysisBlock,
            label = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho25v%s"%i),
            folder = cms.string("CSVIVFv2-NICleanedRho25v%s"%i)    # standard CSVIVFv2 sequence
        )
    )
    tagConfig.append(
        cms.PSet(
            bTagGenericAnalysisBlock,
            label = cms.InputTag("pfDeepCSVJetTagsRho25v%s:probb"%i),
            folder = cms.string("deepCSV_probb-NICleanedRho25v%s"%i)
            )
        )
    # tagConfig.append(
    #     cms.PSet(
    #         bTagGenericAnalysisBlock,
    #         label = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho90v%s"%i),
    #         folder = cms.string("CSVIVFv2-NICleanedRho90v%s"%i)    # standard CSVIVFv2 sequence
    #     )
    # )
for i in range(0, 8):
    tagConfig.append(
        cms.PSet(
            bTagGenericAnalysisBlock,
            label = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTagsCleanedRho9999v%s"%i),
            folder = cms.string("CSVIVFv2-NICleanedRho9999v%s"%i)    # standard CSVIVFv2 sequence
        )
    )
    tagConfig.append(
        cms.PSet(
            bTagGenericAnalysisBlock,
            label = cms.InputTag("pfDeepCSVJetTagsRho9999v%s:probb"%i),
            folder = cms.string("deepCSV_probb-NICleanedRho9999v%s"%i)
            )
        )
    # tagConfig.append(
    #     cms.PSet(
    #         bTagGenericAnalysisBlock,
    #         label = cms.InputTag("pfDeepCSVJetTagsRho90v%s:probb"%i),
    #         folder = cms.string("deepCSV_probb-NICleanedRho90v%s"%i)
    #         )
    #     )

    # tagConfig.append(
    #     cms.PSet(
    #         bTagGenericAnalysisBlock,
    #         label = cms.InputTag("pfDeepFlavourJetTagsRho25v%s:probb"%i),
    #         folder = cms.string("deepJet_probb-NICleanedRho25v%s"%i)
    #         )
    #     )
    # tagConfig.append(
    #     cms.PSet(
    #         bTagGenericAnalysisBlock,
    #         label = cms.InputTag("pfDeepFlavourJetTagsRho90v%s:probb"%i),
    #         folder = cms.string("deepJet_probb-NICleanedRho90v%s"%i)
    #         )
    #     )
    # tagConfig.append(
    #     cms.PSet(
    #         bTagGenericAnalysisBlock,
    #         label = cms.InputTag("pfDeepFlavourJetTagsRho9999v%s:probb"%i),
    #         folder = cms.string("deepJet_probb-NICleanedRho9999v%s"%i)
    #         )
    #     )

"""
end customization
"""

###prints###
#print "Global Tag : ", tag
############

process.load("DQMServices.Components.DQMEnvironment_cfi")
process.load("DQMServices.Core.DQM_cfg")

#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#from RecoBTag.Combined.deepFlavour_cff import *
#from RecoBTag.TensorFlow.pfDeepFlavour_cff import *

#for MC jet flavour
process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")
#process.load("PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi.py")
process.AK4byRef.jets = cms.InputTag("ak4PFJetsCHS")
process.flavourSeq = cms.Sequence(
    process.myPartons *
    process.AK4Flavour
)

process.jetFlavourId = cms.EDProducer("GenJetMatcher",  # cut on deltaR; pick best by deltaR
    src         = cms.InputTag("ak4PFJetsCHS"),      # RECO jets (any View<Jet> is ok)
    matched     = cms.InputTag("ak4GenJets"),        # GEN jets  (must be GenJetCollection)
    mcPdgId     = cms.vint32(),                      # n/a
    mcStatus    = cms.vint32(),                      # n/a
    checkCharge = cms.bool(False),                   # n/a
    maxDeltaR   = cms.double(0.4),                   # Minimum deltaR for the match
    #maxDPtRel   = cms.double(3.0),                  # Minimum deltaPt/Pt for the match (not used in GenJetMatcher)
    resolveAmbiguities    = cms.bool(True),          # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(False),         # False = just match input in order; True = pick lowest deltaR pair first
)

process.patJetPartons = cms.EDProducer('HadronAndPartonSelector',
    src = cms.InputTag("generator"),
    particles = cms.InputTag("genParticles"),
    partonMode = cms.string("Auto"),
    fullChainPhysPartons = cms.bool(True)
)

process.patJetFlavourAssociation = cms.EDProducer("JetFlavourClustering",
    jets = cms.InputTag("ak4PFJetsCHS"),
    bHadrons = cms.InputTag("patJetPartons","bHadrons"),
    cHadrons = cms.InputTag("patJetPartons","cHadrons"),
    partons = cms.InputTag("patJetPartons","physicsPartons"),
    leptons = cms.InputTag("patJetPartons","leptons"),
    jetAlgorithm = cms.string("AntiKt"),
    rParam = cms.double(0.4),
    ghostRescaling = cms.double(1e-18),
    hadronFlavourHasPriority = cms.bool(False)
)


#Validation sequence
process.load("Validation.RecoB.bTagAnalysis_cfi")
process.bTagValidation.jetMCSrc = 'patJetFlavourAssociation'
process.bTagValidation.tagConfig = tagConfig
process.bTagHarvestMC.tagConfig = tagConfig
process.bTagValidation.flavPlots = flavPlots
process.bTagHarvestMC.flavPlots = flavPlots
process.bTagValidation.ptRanges = ptRanges
process.bTagHarvestMC.ptRanges = ptRanges
process.bTagValidation.doPUid = cms.bool(PUid)
process.ak4GenJetsForPUid = cms.EDFilter("GenJetSelector",
                                         src = cms.InputTag("ak4GenJets"),
                                         cut = cms.string('pt > 8.'),
                                         filter = cms.bool(False)
                                         )
process.load("PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi")
process.patJetGenJetMatch.matched = cms.InputTag("ak4GenJetsForPUid")
process.patJetGenJetMatch.maxDeltaR = cms.double(0.25)
process.patJetGenJetMatch.resolveAmbiguities = cms.bool(True)

# load the full reconstraction configuration, to make sure we're getting all needed dependencies
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")

process.load("RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi")
process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:run2_mc")
#process.GlobalTag.globaltag = tag

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

# load b-tagging sequences including NI-rejection

process.load("RecoBTag.SecondaryVertex.bTagNIRejSeqRhoCut25_cff")
process.load("RecoBTag.SecondaryVertex.bTagNIRejSeqRhoCut90_cff")
process.load("RecoBTag.SecondaryVertex.bTagNIRejSeqRhoCut9999_cff")

process.bTagSeq = cms.Sequence(
    process.niRejSeqRhoCut25 *
    process.niRejSeqRhoCut90 *
    process.niRejSeqRhoCut9999
)

#process.deepFlavSeq = cms.Sequence(process.pfDeepFlavourTask)

process.dqmSeq = cms.Sequence( process.ak4GenJetsForPUid * process.patJetGenJetMatch * process.flavourSeq * process.jetFlavourId * process.patJetPartons * process.patJetFlavourAssociation * process.bTagValidation * process.bTagHarvestMC * process.dqmSaver)

#process.plots = cms.Path(process.pfDeepCSV * process.deepFlavSeq * process.bTagSeq * process.dqmSeq)
process.plots = cms.Path(process.bTagSeq * process.dqmSeq)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.customEventContent = cms.PSet(
     outputCommands = cms.untracked.vstring('drop *')
)
process.customEventContent.outputCommands.append('keep *_*_*_validation')

version = '04_TTbar_rho25_rho90_rho999_nirej0to7_aricluspl_deepCSV_beamPipe'
    
process.dqmEnv.subSystemFolder = 'BTAG'
process.dqmSaver.producer = 'DQM'
process.dqmSaver.workflow = '/My/Test/Workflow_%s'%version
process.dqmSaver.convention = 'Offline'
process.dqmSaver.saveByRun = cms.untracked.int32(-1)
process.dqmSaver.saveAtJobEnd =cms.untracked.bool(True) 
process.dqmSaver.forceRunNumber = cms.untracked.int32(1)
process.PoolSource.fileNames = [
#    'file:////nfs/dust/cms/user/eichm/btag/data/2018/RunIISpring18DRPremix_TTToHadronic_TuneCP5_13TeV-powheg-pythia8/328CB6C9-B161-E811-883C-A0369FE2C09C.root'
#    '/store/mc/RunIIAutumn18DRPremix/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/GEN-SIM-RECODEBUG/PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15-v1/80000/FFFD8811-100B-E14D-8AD2-66F035A774F3.root',
#    '/store/mc/RunIIAutumn18DRPremix/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/GEN-SIM-RECODEBUG/PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15-v1/80000/FF2436D6-C3EC-0F4E-8F37-16E953F8F6D8.root'
#    '/store/mc/RunIISpring18DRPremix/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/100X_upgrade2018_realistic_v10-v1/70000/E6845ABF-DE61-E811-A2BF-0663CE00010C.root'
    '/store/mc/RunIISpring18DRPremix/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/100X_upgrade2018_realistic_v10-v1/70000/E0CB7815-2462-E811-9CDA-A0369FE2C146.root'
#    ,'/store/mc/RunIISpring18DRPremix/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/100X_upgrade2018_realistic_v10-v1/70000/D4B1DE5A-B861-E811-B86B-A0369FD0B3B8.root'
#     ,'/store/mc/RunIISpring18DRPremix/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/100X_upgrade2018_realistic_v10-v1/70000/9849E4A5-B161-E811-AC9E-A0369FE2C16A.root'
#     ,'/store/mc/RunIISpring18DRPremix/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/100X_upgrade2018_realistic_v10-v1/70000/54488453-8661-E811-AC57-0CC47A4D99A6.root'
#     ,'/store/mc/RunIISpring18DRPremix/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/100X_upgrade2018_realistic_v10-v1/70000/3CF74428-1362-E811-8CCF-A0369FE2C216.root'
#not existing    ,'/store/mc/RunIISpring18DRPremix/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/100X_upgrade2018_realistic_v10-v1/70000/CE467920-1E62-E811-9F76-A0369FE2C09C.root'
#not existing    ,'/store/mc/RunIISpring18DRPremix/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/100X_upgrade2018_realistic_v10-v1/70000/BCE5B759-8861-E811-B0AC-0CC47A4DEE00.root'
#not existing     ,'/store/mc/RunIISpring18DRPremix/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/100X_upgrade2018_realistic_v10-v1/70000/BAF5F9A2-B761-E811-A09B-A0369FD0B344.root'
#not existing     ,'/store/mc/RunIISpring18DRPremix/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/100X_upgrade2018_realistic_v10-v1/70000/B818B9A2-0E62-E811-BEF6-0CC47A4DECFA.root'
#not existing     ,'/store/mc/RunIISpring18DRPremix/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/100X_upgrade2018_realistic_v10-v1/70000/9096A7A4-1D62-E811-BC13-A0369FE2C19A.root'
#not existing     ,'/store/mc/RunIISpring18DRPremix/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/100X_upgrade2018_realistic_v10-v1/70000/8A77C800-DF61-E811-923C-A0369FE2C0D0.root'
#not existing     ,'/store/mc/RunIISpring18DRPremix/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/100X_upgrade2018_realistic_v10-v1/70000/328CB6C9-B161-E811-883C-A0369FE2C09C.root'
#    'file:////nfs/dust/cms/user/eichm/btag/data/TTToHadronic_RunIISpring18DRPremix_deepJet_E6845ABF-DE61-E811-A2BF-0663CE00010C.root'
]





outFile = open("tmpConfig_validation.py","w")
outFile.write(process.dumpPython())
outFile.close()

#  LocalWords:  PFJetsMCFlavourInfos
