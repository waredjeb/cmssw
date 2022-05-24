import FWCore.ParameterSet.Config as cms

from RecoBTag.SecondaryVertex.nuclearInteractionIdentifier_cfi import *
from RecoBTag.SecondaryVertex.vertexAndTracksCleaner_cfi import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
from RecoBTag.SecondaryVertex.trackRefitter_cfi import *



#===========================================
# NI rejection and tracks cleaning procedure
#
#===========================================


processName = "DATA"
process = cms.Process(processName)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source ("PoolSource",
    fileNames=cms.untracked.vstring(
#        'file:////nfs/dust/cms/user/eichm/btag/data/2018/RunIISpring18DRPremix_TTToHadronic_TuneCP5_13TeV-powheg-pythia8/328CB6C9-B161-E811-883C-A0369FE2C09C.root'
#        'file:////pnfs/desy.de/cms/tier2/store/data/Run2017B/SingleMuon/MINIAOD/PromptReco-v2/000/299/000/00000/7647626E-2D6A-E711-829C-02163E01A606.root'
#        'file:////pnfs/desy.de/cms/tier2/store/data/Run2016C/SingleMuon/AOD/23Sep2016-v1/80002/8CEDE132-C78A-E611-962D-002590D9D980.root'
        'file:////pnfs/desy.de/cms/tier2/store/data/Run2017B/SingleMuon/AOD/17Nov2017-v1/70008/0689C2B4-37D8-E711-ACAF-02163E01A208.root'
    )
)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v6'
#process.GlobalTag = GlobalTag(process.GlobalTag, "auto:run2_mc")

# JSON filter
import FWCore.PythonUtilities.LumiList as LumiList
jsonName = "Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON" #golden json
process.source.lumisToProcess = LumiList.LumiList(filename = 'data/'+jsonName+'.txt').getVLuminosityBlockRange()
print "JSON file loaded: ", jsonName

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
    primaryVertices  = cms.InputTag("offlineSlimmedPrimaryVertices"),
    secondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
    selection = cms.PSet(
        nuclearInteractionCandIdentifier.selection,
        position = cms.vdouble(2.17, 2.25, 2.694, 3.481, 6.597, 7.278, 10.700, 11.344, 15.797, 16.430)
    )
)
#process.vertexRefitted.secondaryVertices = cms.InputTag("nuclearInteractionIdentifier0")
process.vertexRefitted0 = cms.EDProducer("NITrackReFitter",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    secondaryVertices = cms.InputTag("nuclearInteractionIdentifier0")
)

process.nuclearInteractionIdentifierAfterRefit = cms.EDProducer("NuclearInteractionCandidateIdentifier",
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
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
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
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
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    secondaryVertices = cms.InputTag("candidateVertexMerger"),
    sigCut = cms.double(5),
    trackMinLayers = cms.int32(4),
    trackMinPixels = cms.int32(1),
    trackMinPt = cms.double(0.4),
    tracks = cms.InputTag("vertexAndTracksCandCleaned0")
)

process.slimmedSecondaryVertices0 = cms.EDProducer("CandidateVertexMerger",
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



process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string("/nfs/dust/cms/user/eichm/btag/ntuple/vertex_QCD_170to300_3.root"),
#    fileName = cms.untracked.string("/store/user/meich/NI_identification/vertex_QCD_default_test.root"),
#    fileName = cms.untracked.string("/nfs/dust/cms/user/eichm/btag/ntuple/DATA2017.root"),
     fileName = cms.untracked.string("DATA2017.root"),
    outputCommands = cms.untracked.vstring("drop *",
                        "keep *_*_*_RECODEBUG",
                        "keep *_*_*_SIM",
                        "keep *_inclusiveSecondaryVertices_*_*",
                        "keep *_slimmedSecondaryVertices_*_*",
                        "keep *_offlineSlimmedPrimaryVertices_*_*",
                        "keep *_offlinePrimaryVertices_*_*",
                        "keep *_inclusiveCandidateSecondaryVertices_*_*",
                        "keep *_particleFlowDisplacedVertex_*_*",
                        "keep *_*Track*_*_*",
                        "keep *_*Vertex*_*_*",
                        "keep *_*Vertices*_*_*",
                        "keep *_mix_MergedTrackTruth_*")
)

process.p = cms.Path(process.inclusiveCandidateVertexing)
# * process.nuclearInteractionIdentifier0 )
#* process.vertexRefitted0 * process.nuclearInteractionIdentifierAfterRefit )
#* process.vertexAndTracksCandCleaned0 * process.inclusiveCandidateVertexFinder0 * process.candidateVertexMerger0 *  process.candidateVertexArbitrator0 * process.inclusiveSecondaryVerticesCleaned0)

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
