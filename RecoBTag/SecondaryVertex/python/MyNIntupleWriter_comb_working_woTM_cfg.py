import FWCore.ParameterSet.Config as cms

from RecoBTag.SecondaryVertex.nuclearInteractionIdentifier_cfi import *
from RecoBTag.SecondaryVertex.vertexAndTracksCleaner_cfi import *
#from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
from RecoVertex.AdaptiveVertexFinder.inclusiveCandidateVertexFinder_cfi import *
from RecoBTag.SecondaryVertex.trackRefitter_cfi import *



#===========================================
# NI rejection and tracks cleaning procedure
#
#===========================================


processName = "comb"
process = cms.Process(processName)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source ("PoolSource",
    fileNames=cms.untracked.vstring(
        'file:////nfs/dust/cms/user/eichm/btag/data/2018/RunIISpring18DRPremix_TTToHadronic_TuneCP5_13TeV-powheg-pythia8/328CB6C9-B161-E811-883C-A0369FE2C09C.root'
#'file:////nfs/dust/cms/user/eichm/btag/ntuple/NITMTTbar18.root'
    )
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = '94X_mc2017_realistic_v10'
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
process.vertexAndTracksCandCleaned.veto = cms.InputTag('nuclearInteractionCandIdentifier0')
#process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveCandidateVertexFinder_cfi")
process.load("RecoBTag.SecondaryVertex.trackRefitter_cfi")


process.nuclearInteractionIdentifier0 = cms.EDProducer("NuclearInteractionCandidateIdentifier",
    primaryVertices  = cms.InputTag("offlinePrimaryVertices"),
    secondaryVertices = cms.InputTag("inclusiveCandidateSecondaryVertices"),
    selection = cms.PSet(
        nuclearInteractionCandIdentifier.selection,
        position = cms.vdouble(2.17, 2.25, 2.694, 3.481, 6.597, 7.278, 10.700, 11.344, 15.797, 16.430)
    )
)
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
    maximumLongitudinalImpactParameter = cms.double(99999.9),
    maximumTimeSignificance = cms.double(99999.9),
    minHits = cms.uint32(1),
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
    secondaryVertices = cms.InputTag("candidateVertexMerger0"),
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

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("SimTracker.TrackHistory.VertexClassifier_cff")
process.vertexClassifier.vertexProducer = cms.untracked.InputTag("offlinePrimaryVertices")
process.vertexHistory.vertexProducer = cms.untracked.InputTag("offlinePrimaryVertices")

process.ak4JetTracksAssociatorAtVertexPF = cms.EDProducer("JetTracksAssociatorAtVertex",
    tracks = cms.InputTag("generalTracks"),
    coneSize = cms.double(0.4),
    useAssigned = cms.bool(True),
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    jets = cms.InputTag("ak4PFJetsCHS")
)


process.load("SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi")
process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")

process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)

process.load("RecoBTag.ImpactParameter.pfImpactParameterTagInfos_cfi")
process.pfImpactParameterTagInfos.primaryVertex = cms.InputTag("offlinePrimaryVertices")
process.pfImpactParameterTagInfos.candidates = ("vertexAndTracksCandCleaned0")

# SecondaryVertexProducer using IVF Vertex collection
process.load("RecoBTag.SecondaryVertex.pfInclusiveSecondaryVertexFinderTagInfos_cfi")
process.pfInclusiveSecondaryVertexFinderTagInfos.extSVCollection = cms.InputTag("inclusiveSecondaryVerticesCleaned0")
process.pfInclusiveSecondaryVertexFinderTagInfos.vertexCuts.distVal2dMax = cms.double(99999.9)
process.pfInclusiveSecondaryVertexFinderTagInfos.vertexCuts.distSig2dMin = cms.double(-99999.9)
process.vertexCutsBlock.vertexCuts.distVal2dMax = cms.double(99999.9)
process.vertexCutsBlock.vertexCuts.distSig2dMin = cms.double(-99999.9)
process.pfInclusiveSecondaryVertexFinderTagInfos.useExternalSV = cms.bool(True)
process.pfInclusiveSecondaryVertexFinderTagInfos.useSVClustering = cms.bool(True)
process.pfInclusiveSecondaryVertexFinderTagInfos.jetAlgorithm = cms.string("AntiKt")
process.pfInclusiveSecondaryVertexFinderTagInfos.rParam = cms.double(0.4)

process.load("SimTracker.TrackHistory.SecondaryVertexTagInfoProxy_cff")
process.svTagInfoProxy.svTagInfoProducer = cms.untracked.InputTag("pfInclusiveSecondaryVertexFinderTagInfos")

process.load("PhysicsTools.UtilAlgos.TFileService_cfi")

process.svTagInfoValidation = cms.EDAnalyzer("SVTagInfoValidationAnalyzer",
    svTagInfoProducer = cms.untracked.InputTag("pfInclusiveSecondaryVertexFinderTagInfos"),
    bestMatchByMaxValue = cms.untracked.bool(True),
    trackingTruth = cms.untracked.InputTag('mix','MergedTrackTruth'),
    vertexAssociator = cms.untracked.InputTag('vertexAssociatorByTracksByHits'),
    vertexProducer = cms.untracked.InputTag('svTagInfoProxy'),
    enableRecoToSim = cms.untracked.bool(True),
    enableSimToReco = cms.untracked.bool(False),
    hepMC = cms.untracked.InputTag("generatorSmeared"),
    longLivedDecayLength = cms.untracked.double(1e-14),
    vertexClusteringDistance = cms.untracked.double(0.003)
)

process.svTagInfoValidationNImatch = cms.EDAnalyzer("SVTagInfoValidationNImatchAnalyzer",
    svTagInfoProducer = cms.untracked.InputTag("pfInclusiveSecondaryVertexFinderTagInfos"),
    secondaryVertices = cms.untracked.InputTag("inclusiveSecondaryVerticesCleaned0"),
#    secondaryVertices = cms.untracked.InputTag("inclusiveCandidateSecondaryVertices"),
    bestMatchByMaxValue = cms.untracked.bool(True),
    trackingTruth = cms.untracked.InputTag('mix','MergedTrackTruth'),
    vertexAssociator = cms.untracked.InputTag('vertexAssociatorByTracksByHits'),
    vertexProducer = cms.untracked.InputTag('svTagInfoProxy'),
    enableRecoToSim = cms.untracked.bool(True),
    enableSimToReco = cms.untracked.bool(False),
    hepMC = cms.untracked.InputTag("generatorSmeared"),
    longLivedDecayLength = cms.untracked.double(1e-14),
    vertexClusteringDistance = cms.untracked.double(0.003)
)

process.trackingParticleRecoTrackAsssociationByHits.associator = cms.InputTag("quickTrackAssociatorByHits")

process.vertexHistoryAnalyzer = cms.EDAnalyzer("VertexHistoryAnalyzer",
    process.vertexClassifier
)






process.out = cms.OutputModule("PoolOutputModule",
#     fileName = cms.untracked.string("/nfs/dust/cms/user/eichm/btag/ntuple/NITTbar18_test.root"),
     fileName = cms.untracked.string("/nfs/dust/cms/user/eichm/btag/ntuple/NIcombInfo_test.root"),
    outputCommands = cms.untracked.vstring("keep *",
#                        "keep *_*_*_*")
                        "keep *_*_*_IdAndTM",
                        "keep *_*_*_SIM",
                        "keep *_inclusiveSecondaryVertices_*_*",
                        "keep *_offlinePrimaryVertices_*_*",
                        "keep *_offlineBeamSpot_*_*",
                        "keep *GenParticle*_*_*_*",
                        "keep *Jet*_*_*_*",
                        "keep *Candidate*_*_*_*",
                        "keep *_generalTracks_*_*",
                        "keep *_mix_MergedTrackTruth_*")
)


process.nuclearIdentification = cms.Sequence(process.inclusiveCandidateVertexFinder * process.nuclearInteractionIdentifier0 * process.vertexRefitted0 * process.nuclearInteractionIdentifierAfterRefit * process.vertexAndTracksCandCleaned0 * process.inclusiveCandidateVertexFinder0 * process.candidateVertexMerger0 *  process.candidateVertexArbitrator0 * process.inclusiveSecondaryVerticesCleaned0 *
 process.tpClusterProducer * process.quickTrackAssociatorByHits* process.trackingParticleRecoTrackAsssociationByHits * process.ak4JetTracksAssociatorAtVertexPF * process.vertexAssociatorByTracksByHits * process.pfImpactParameterTagInfos * process.pfInclusiveSecondaryVertexFinderTagInfos )
#* process.svTagInfoProxy * process.svTagInfoValidationNImatch)

process.p = cms.Path(process.nuclearIdentification )

process.e = cms.EndPath(process.out)

outFile = open("tmpConfig_comb.py","w")
outFile.write(process.dumpPython())
outFile.close()
