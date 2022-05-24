import FWCore.ParameterSet.Config as cms

from RecoBTag.SecondaryVertex.nuclearInteractionIdentifier_cfi import *
from RecoBTag.SecondaryVertex.vertexAndTracksCleaner_cfi import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
from RecoBTag.SecondaryVertex.trackRefitter_cfi import *



#===========================================
# NI rejection and tracks cleaning procedure
#
#===========================================


processName = "RECODEBUG"
process = cms.Process(processName)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source ("PoolSource",
    fileNames=cms.untracked.vstring(
        'file:////nfs/dust/cms/user/eichm/btag/data/2018/RunIISpring18DRPremix_TTToHadronic_TuneCP5_13TeV-powheg-pythia8/328CB6C9-B161-E811-883C-A0369FE2C09C.root'
#        'file:////nfs/dust/cms/user/eichm/btag/data/bgd/84B70053-DA98-E711-8525-00259029E87C.root'
#         'file:////nfs/dust/cms/user/eichm/btag/data/2016/00104622-EA3C-E611-807A-002590200828.root'
#        'file:////nfs/dust/cms/user/eichm/btag/data/2017/RunIIFall17DRPremix_QCD_Pt_50to80_TuneCP5_13TeV_pythia8/4299C2F1-4FDC-E711-938D-A4BF0112BC52.root'
# #        '/store/mc/RunIIFall17DRPremix/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/GEN-SIM-RECODEBUG/94X_mc2017_realistic_v10-v1/40000/008D9025-4CDC-E711-A398-0242AC130002.root'
#        'file:////nfs/dust/cms/user/eichm/btag/data/2017/RunIIFall17DRPremix_QCD_Pt_80to120_TuneCP5_13TeV_pythia8/72CD356F-4BDC-E711-885D-0242AC130002.root'
#        '/store/mc/RunIIFall17DRPremix/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/GEN-SIM-RECODEBUG/94X_mc2017_realistic_v10-v1/40000/7A8DF3B6-E2DB-E711-BEEB-0242AC130002.root'
#         '/store/mc/RunIIFall17DRPremix/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/GEN-SIM-RECODEBUG/94X_mc2017_realistic_v10-v1/40000/BC680365-52DC-E711-B825-FA163E04137F.root'
# #        '/store/mc/RunIIFall17DRPremix/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/GEN-SIM-RECODEBUG/94X_mc2017_realistic_v10-v1/40000/CA84E2A4-4FDC-E711-A204-FA163EC5DEF2.root'
    )
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = '94X_mc2017_realistic_forppRef5TeV'
#process.GlobalTag.globaltag = '94X_mc2017_realistic_v10'
#process.GlobalTag = GlobalTag(process.GlobalTag, "auto:run2_mc")

process.load("TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi")
process.load("TrackingTools.TrackFitters.TrackFitters_cff")
process.load("TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi")
process.load("TrackingTools.KalmanUpdators.KFUpdatorESProducer_cfi")


process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")

process.load("RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi")
process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder') 
)

process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff")

process.load("RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff")

process.load("RecoBTag.SecondaryVertex.vertexAndTracksCleaner_cfi")
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
process.load("RecoBTag.SecondaryVertex.trackRefitter_cfi")

#process.load("RecoBTag.SecondaryVertex.nuclearInteractionIdentifier_cfi")
# process.nuclearInteractionCandIdentifier.selection.position = cms.vdouble(2.17, 2.25, 2.694, 3.481, 6.597, 7.278, 10.700, 11.344, 15.797, 16.430)

process.nuclearInteractionIdentifier0 = cms.EDProducer("NuclearInteractionCandidateIdentifier",
    primaryVertices  = cms.InputTag("offlinePrimaryVertices"),
    secondaryVertices = cms.InputTag("inclusiveCandidateSecondaryVertices"),
    selection = cms.PSet(
        nuclearInteractionCandIdentifier.selection,
        position = cms.vdouble(2.17, 2.25, 2.694, 3.481, 6.597, 7.278, 10.700, 11.344, 15.797, 16.430)
    )
)
#process.vertexRefitted.secondaryVertices = cms.InputTag("nuclearInteractionIdentifier0")
process.vertexRefitted0 = cms.EDProducer("NITrackReFitter",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    secondaryVertices = cms.InputTag("nuclearInteractionIdentifier0")
)

process.nuclearInteractionIdentifierAfterRefit = cms.EDProducer("NuclearInteractionCandidateIdentifier",
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    secondaryVertices = cms.InputTag("vertexRefitted0"),
    selection = cms.PSet(
    position = cms.vdouble(2.17, 2.25, 2.694, 3.481, 6.597, 
        7.278, 10.7, 11.344, 15.797, 16.43)
    )
)

#process.vertexAndTracksCandCleaned.veto = cms.InputTag("nuclearInteractionIdentifierAfterRefit")
#process.inclusiveCandidateVertexFinder.tracks = cms.InputTag("vertexAndTracksCandCleaned")
#process.candidateVertexMerger.secondaryVertices = cms.InputTag("inclusiveCandidateVertexFinder")
#process.candidateVertexArbitrator.tracks = cms.InputTag("vertexAndTracksCandCleaned")
#process.candidateVertexArbitrator.secondaryVertices = cms.InputTag("candidateVertexMerger")

process.inclusiveCandidateVertexFinder0 = cms.EDProducer("InclusiveCandidateVertexFinder",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clusterizer = cms.PSet(
        clusterMaxDistance = cms.double(0.05),
        clusterMaxSignificance = cms.double(4.5),
        clusterMinAngleCosine = cms.double(0.5),
        distanceRatio = cms.double(20),
        maxTimeSignificance = cms.double(3.5),
        seedMax3DIPSignificance = cms.double(9999),
        seedMax3DIPValue = cms.double(9999),
        seedMin3DIPSignificance = cms.double(1.2),
        seedMin3DIPValue = cms.double(0.005)
    ),
    fitterRatio = cms.double(0.25),
    fitterSigmacut = cms.double(3),
    fitterTini = cms.double(256),
    maxNTracks = cms.uint32(30),
    maximumLongitudinalImpactParameter = cms.double(0.3),
    maximumTimeSignificance = cms.double(3),
    minHits = cms.uint32(0),
    minPt = cms.double(0.8),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    tracks = cms.InputTag("vertexAndTracksCandCleaned0"),
    useDirectVertexFitter = cms.bool(True),
    useVertexReco = cms.bool(True),
    vertexMinAngleCosine = cms.double(0.95),
    vertexMinDLen2DSig = cms.double(2.5),
    vertexMinDLenSig = cms.double(0.5),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        primcut = cms.double(1),
        seccut = cms.double(3),
        smoothing = cms.bool(True)
    )
)

process.vertexAndTracksCandCleaned0 = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("particleFlow"),
    veto = cms.InputTag("nuclearInteractionIdentifierAfterRefit")
)

process.candidateVertexArbitrator0 = cms.EDProducer("CandidateVertexArbitrator",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    dLenFraction = cms.double(0.333),
    dRCut = cms.double(0.4),
    distCut = cms.double(0.04),
    fitterRatio = cms.double(0.25),
    fitterSigmacut = cms.double(3),
    fitterTini = cms.double(256),
    maxTimeSignificance = cms.double(3.5),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    secondaryVertices = cms.InputTag("candidateVertexMerger"),
    sigCut = cms.double(5),
    trackMinLayers = cms.int32(4),
    trackMinPixels = cms.int32(1),
    trackMinPt = cms.double(0.4),
    tracks = cms.InputTag("vertexAndTracksCandCleaned0")
)

process.inclusiveCandidateSecondaryVertices0 = cms.EDProducer("CandidateVertexMerger",
    maxFraction = cms.double(0.2),
    minSignificance = cms.double(10.0),
    secondaryVertices = cms.InputTag("candidateVertexArbitrator0")
)


process.candidateVertexMerger0 = cms.EDProducer("CandidateVertexMerger",
    maxFraction = cms.double(0.7),
    minSignificance = cms.double(2),
    secondaryVertices = cms.InputTag("inclusiveCandidateVertexFinder0")
)

process.inclusiveSecondaryVerticesCleaned0 = cms.EDProducer("CandidateVertexMerger",
    maxFraction = cms.double(0.2),
    minSignificance = cms.double(10.),
    secondaryVertices = cms.InputTag("candidateVertexArbitrator0")
)



# process.out = cms.OutputModule("PoolOutputModule",
#     fileName = cms.untracked.string("/nfs/dust/cms/user/eichm/btag/ntuple/NI_TTbar18.root")
# )

process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string("/nfs/dust/cms/user/eichm/btag/ntuple/vertex_QCD_170to300_3.root"),
#    fileName = cms.untracked.string("/store/user/meich/NI_identification/vertex_QCD_default_test.root"),
    fileName = cms.untracked.string("/nfs/dust/cms/user/eichm/btag/ntuple/NITTbar18_test.root"),
    outputCommands = cms.untracked.vstring("keep *",
                        "keep *_*_*_*")
                        # "keep *_*_*_RECODEBUG",
                        # "keep *_*_*_SIM",
                        # "keep *_inclusiveSecondaryVertices_*_*",
                        # "keep *_offlinePrimaryVertices_*_*",
                        # "keep *GenParticle*_*_*_*",
                        # "keep *Jet*_*_*_*",
                        # "keep *Candidate*_*_*_*",
                        # "keep *_generalTracks_*_*",
                        # "keep *_mix_MergedTrackTruth_*")
)

process.p = cms.Path(process.inclusiveCandidateVertexing * process.nuclearInteractionIdentifier0 * process.vertexRefitted0 * process.nuclearInteractionIdentifierAfterRefit * process.vertexAndTracksCandCleaned0 * process.inclusiveCandidateVertexFinder0 * process.candidateVertexMerger0 *  process.candidateVertexArbitrator0 * process.inclusiveSecondaryVerticesCleaned0)

#process.p = cms.Path(process.nuclearInteractionIdentifier0 * process.vertexRefitted * process.nuclearInteractionIdentifierAfterRefit * process.vertexAndTracksCandCleaned * process.inclusiveCandidateVertexFinder * process.candidateVertexMerger * process.candidateVertexArbitrator * process.inclusiveSecondaryVerticesCleaned)
process.e = cms.EndPath(process.out)

# all NI rejection sequences

# nuclearInteractionsRemoved0 = cms.Sequence(
# 	inclusiveCandidateVertexing *
# 	nuclearInteractionCandIdentifier *
# 	vertexRefitted *
# 	nuclearInteractionIdentifierAfterRefit *
# 	vertexAndTracksCleaned *
# 	inclusiveVertexFinderCleaned *
# 	vertexMergerCleaned *
# 	trackVertexArbitratorCleaned *
# 	inclusiveSecondaryVerticesCleaned
# )

outFile = open("tmpConfig_NICand.py","w")
outFile.write(process.dumpPython())
outFile.close()
