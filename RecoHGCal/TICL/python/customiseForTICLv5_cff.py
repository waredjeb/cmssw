import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.ticlDumper_cfi import ticlDumper
from RecoHGCal.Configuration.RecoHGCal_EventContent_cff import customiseForTICLv5EventContent
from RecoHGCal.TICL.ticlDumper_cfi import ticlDumper
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterFromCPsAssociationPR as _tracksterSimTracksterFromCPsAssociationPR
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterAssociationPR  as _tracksterSimTracksterAssociationPR

def customiseTICLv5FromReco(process, enableDumper = False):
    # TensorFlow ESSource
    process.TFESSource = cms.Task(process.trackdnn_source)

    # Reconstruction
    process.hgcalLayerClustersTask = cms.Task(process.hgcalLayerClustersEE,
                                              process.hgcalLayerClustersHSi,
                                              process.hgcalLayerClustersHSci,
                                              process.hgcalMergeLayerClusters)

    process.ticlIterationsTask = cms.Task(
        process.ticlCLUE3DHighStepTask,
        process.ticlPassthroughStepTask,
        process.ticlTracksterLinksTask
    )

    process.mergeTICLTask = cms.Task()

    process.iterTICLTask = cms.Path(process.hgcalLayerClustersTask,
                            process.TFESSource,
                            process.ticlLayerTileTask,
 #                           process.mtdSoATask,
                            process.mergeTICLTask,
                            process.ticlIterationsTask,
                            process.ticlCandidateTask,
                            process.ticlPFTask)

    process.tracksterSimTracksterFromCPsAssociationPRHigh = _tracksterSimTracksterFromCPsAssociationPR.clone(
        label_tst = cms.InputTag("ticlTrackstersCLUE3DHigh")
        )
    process.tracksterSimTracksterAssociationPRHigh = _tracksterSimTracksterAssociationPR.clone(
        label_tst = cms.InputTag("ticlTrackstersCLUE3DHigh")
        )

    process.hgcalAssociators = cms.Task(process.recHitMapProducer, process.lcAssocByEnergyScoreProducer, process.layerClusterCaloParticleAssociationProducer,
                            process.scAssocByEnergyScoreProducer, process.layerClusterSimClusterAssociationProducer,
                            # FP 07/2024 new associators:
                            layerClusterToCLUE3DTracksterAssociation, layerClusterToTracksterMergeAssociation,
                            layerClusterToSimTracksterAssociation,
                            hitToTrackstersAssociationLinking, hitToTrackstersAssociationPR,
                            hitToSimTracksterAssociation, hitToSimTracksterFromCPsAssociation,
                            tracksterSimTracksterAssociationByHitsLinking, tracksterSimTracksterAssociationByHitsPR,
                            tracksterSimTracksterFromCPsAssociationLinking, tracksterSimTracksterAssociationLinking, tracksterSimTracksterFromCPsAssociationPR, tracksterSimTracksterAssociationPR,
                            hitToSimClusterCaloParticleAssociator
                            )

    if(enableDumper):
        process.ticlDumper = ticlDumper.clone(
            saveLCs=True,
            saveCLUE3DTracksters=True,
            saveTrackstersMerged=True,
            saveSimTrackstersSC=True,
            saveSimTrackstersCP=True,
            saveTICLCandidate=True,
            saveSimTICLCandidate=True,
            saveTracks=True,
            saveAssociations=True,
            trackstersclue3d = cms.InputTag('ticlTrackstersCLUE3DHigh'),
            ticlcandidates = cms.InputTag("ticlCandidate"),
            trackstersmerged = cms.InputTag("ticlCandidate"),
            trackstersInCand = cms.InputTag("ticlCandidate")
        )
        process.TFileService = cms.Service("TFileService",
                                           fileName=cms.string("histo.root")
                                           )

        process.FEVTDEBUGHLToutput_step = cms.EndPath(process.ticlDumper)

    process.TICL_Validator = cms.Task(process.hgcalValidator)
    process.TICL_Validation = cms.Path(process.ticlSimTrackstersTask, process.hgcalAssociators, process.TICL_Validator)

    # Schedule definition
    process.schedule = cms.Schedule(process.iterTICLTask,
                                    process.TICL_Validation,
                                    process.FEVTDEBUGHLToutput_step)
    process = customiseForTICLv5EventContent(process)

    return process
