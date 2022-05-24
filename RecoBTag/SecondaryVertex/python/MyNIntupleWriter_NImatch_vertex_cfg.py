import FWCore.ParameterSet.Config as cms
import os
import copy

import SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi as tabh

processName = "TruthMatching"
process = cms.Process(processName)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

## Messagge logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    allowUnscheduled = cms.untracked.bool(True),
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("SimTracker.TrackHistory.VertexClassifier_cff")
process.vertexClassifier.vertexProducer = cms.untracked.InputTag("vertexRefitted0")
process.vertexHistory.vertexProducer = cms.untracked.InputTag("vertexRefitted0") 

process.ak4JetTracksAssociatorAtVertexPF = cms.EDProducer("JetTracksAssociatorAtVertex",
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.4),
    useAssigned = cms.bool(True),
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    jets = cms.InputTag("ak4PFJetsCHS")
)


process.load("SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi")
process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:run2_mc")

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")


process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)

process.load("RecoBTag.ImpactParameter.impactParameterTagInfos_cfi")
process.impactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
process.impactParameterTagInfos.maximumLongitudinalImpactParameter = cms.double(99999.9)
process.impactParameterTagInfos.maximumChiSquared  = cms.double(99999.9)
process.impactParameterTagInfos.maximumTransverseImpactParameter = cms.double(99999.9)

# SecondaryVertexProducer using IVF Vertex collection
process.load("RecoBTag.SecondaryVertex.secondaryVertexTagInfos_cfi")
#process.secondaryVertexTagInfos.extSVCollection = cms.InputTag("vertexRefitted0")
# process.secondaryVertexTagInfos.vertexCuts.distVal2dMax = cms.double(2.5)
# process.secondaryVertexTagInfos.vertexCuts.distSig2dMin = cms.double(3.0)
# process.vertexCutsBlock.vertexCuts.distVal2dMax = cms.double(2.5)
# process.vertexCutsBlock.vertexCuts.distSig2dMin = cms.double(3.0)
process.secondaryVertexTagInfos.vertexCuts.distVal2dMax = cms.double(99999.9)
process.secondaryVertexTagInfos.vertexCuts.distSig2dMin = cms.double(-99999.9)
process.vertexCutsBlock.vertexCuts.distVal2dMax = cms.double(99999.9)
process.vertexCutsBlock.vertexCuts.distSig2dMin = cms.double(-99999.9)
process.secondaryVertexTagInfos.useExternalSV = cms.bool(False)
process.secondaryVertexTagInfos.useSVClustering = cms.bool(True)
process.secondaryVertexTagInfos.jetAlgorithm = cms.string("AntiKt")
process.secondaryVertexTagInfos.rParam = cms.double(0.4)

process.load("SimTracker.TrackHistory.SecondaryVertexTagInfoProxy_cff")
process.svTagInfoProxy.svTagInfoProducer = cms.untracked.InputTag("secondaryVertexTagInfos")

process.load("PhysicsTools.UtilAlgos.TFileService_cfi")

process.svTagInfoValidation = cms.EDAnalyzer("SVTagInfoValidationAnalyzer",
    svTagInfoProducer = cms.untracked.InputTag("secondaryVertexTagInfos"),
    bestMatchByMaxValue = cms.untracked.bool(True),
    trackingTruth = cms.untracked.InputTag('mix','MergedTrackTruth'),
    vertexAssociator = cms.untracked.InputTag('vertexAssociatorByTracksByHits'),
    vertexProducer = cms.untracked.InputTag('svTagInfoProxy'),
    enableRecoToSim = cms.untracked.bool(True),
    enableSimToReco = cms.untracked.bool(False),
    hepMC = cms.untracked.InputTag("generatorSmeared"),
    longLivedDecayLength = cms.untracked.double(1e-14),
    vertexClusteringDistance = cms.untracked.double(0.003),
    jets = cms.untracked.InputTag("ak4PFJetsCHS")
)

process.trackingParticleRecoTrackAsssociationByHits.associator = cms.InputTag("quickTrackAssociatorByHits")

process.vertexHistoryAnalyzer = cms.EDAnalyzer("VertexHistoryAnalyzer",
    process.vertexClassifier
)


# source for TTbar Spring2018 sample
process.source = cms.Source ("PoolSource",
    fileNames=cms.untracked.vstring(
        '/store/mc/RunIIFall17DRPremix/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/GEN-SIM-RECODEBUG/94X_mc2017_realistic_v10-v1/40000/CA84E2A4-4FDC-E711-A204-FA163EC5DEF2.root'
#        ' /store/mc/RunIIFall17DRPremix/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/GEN-SIM-RECODEBUG/94X_mc2017_realistic_v10-v1/50000/125C3D2F-68E2-E711-AB6E-E0071B695B81.root'
#        '/store/mc/RunIIFall17DRPremix/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/GEN-SIM-RECODEBUG/PU2017_94X_mc2017_realistic_v11_ext1-v2/20000/00B8E806-4C5E-E811-9D01-5065F381E182.root'
#        '/store/mc/RunIIAutumn18DRPremix/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/GEN-SIM-RECODEBUG/PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15-v1/80000/FFFD8811-100B-E14D-8AD2-66F035A774F3.root'
#        'file:////nfs/dust/cms/user/eichm/btag/data/2018/RunIISpring18DRPremix_TTToHadronic_TuneCP5_13TeV-powheg-pythia8/328CB6C9-B161-E811-883C-A0369FE2C09C.root'
#        'file:////nfs/dust/cms/user/eichm/btag/ntuple/NITTbar18_test.root'
    )
)

process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string("/nfs/dust/cms/user/eichm/btag/ntuple/vertex_TTbar18_test.root"),
#    fileName = cms.untracked.string("/nfs/dust/cms/user/eichm/btag/ntuple/vertex_QCD_80to120_v11_test.root"),
    fileName = cms.untracked.string("/nfs/dust/cms/user/eichm/btag/ntuple/vertex_TTbar18_withoutCut.root"),
    outputCommands = cms.untracked.vstring("keep *",
                        "keep *_*_*_TruthMatching",
                        "keep *_inclusiveSecondaryVertices_*_*",
                        "keep *_offlinePrimaryVertices_*_*",
                        "keep *_*_*_SIM",
                        "keep *GenParticle*_*_*_*",
                        "keep *PFJet*_*_*_*",
                        "keep *_*PFJet*_*_*",
                        "keep *Candidate*_*_*_*",
                        "keep *_generalTracks_*_*",
                        "keep *_mix_MergedTrackTruth_*")
)

process.vertexAssociatorSequence = cms.Sequence( process.tpClusterProducer * process.quickTrackAssociatorByHits* process.trackingParticleRecoTrackAsssociationByHits * process.ak4JetTracksAssociatorAtVertexPF * process.vertexAssociatorByTracksByHits * process.impactParameterTagInfos * process.secondaryVertexTagInfos )

process.p = cms.Path( process.vertexAssociatorSequence  * process.svTagInfoProxy * process.svTagInfoValidation)

process.e = cms.EndPath(process.out)

outFile = open("tmpConfig_vertex.py","w")
outFile.write(process.dumpPython())
outFile.close()
