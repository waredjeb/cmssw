# Reconstruction
from RecoHGCal.TICL.iterativeTICL_cff import *
from RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cff import hgcalLayerClustersEE, hgcalLayerClustersHSi, hgcalLayerClustersHSci
from RecoLocalCalo.HGCalRecProducers.hgcalMergeLayerClusters_cfi import hgcalMergeLayerClusters
from RecoHGCal.TICL.ticlDumper_cfi import ticlDumper
from RecoHGCal.TICL.ticlGraphAnalyzer_cfi import ticlGraphAnalyzer
# Validation
from Validation.HGCalValidation.HGVHistoProducerAlgoBlock_cfi import *
from Validation.HGCalValidation.HGCalValidator_cfi import *
from RecoLocalCalo.HGCalRecProducers.hgcalRecHitMapProducer_cfi import hgcalRecHitMapProducer

# Load DNN ESSource
from RecoTracker.IterativeTracking.iterativeTk_cff import trackdnn_source
from RecoHGCal.TICL.ticlGraphProducer_cfi import ticlGraphProducer as _ticlGraphProducer
from RecoHGCal.TICL.SimTracksters_cff import *
from RecoHGCal.TICL.simpleValidation_cfi import *
# Automatic addition of the customisation function from RecoHGCal.Configuration.RecoHGCal_EventContent_cff
from RecoHGCal.Configuration.RecoHGCal_EventContent_cff import customiseHGCalOnlyEventContent
from SimCalorimetry.HGCalAssociatorProducers.simTracksterAssociatorByEnergyScore_cfi import simTracksterAssociatorByEnergyScore as simTsAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterAssociationLinking, tracksterSimTracksterAssociationPR, tracksterSimTracksterAssociationLinkingbyCLUE3D, tracksterSimTracksterAssociationPRbyCLUE3D, tracksterSimTracksterAssociationLinkingPU, tracksterSimTracksterAssociationPRPU,tracksterSimTracksterAssociationLinkingbyCLUE3DPU,tracksterSimTracksterAssociationPRbyCLUE3DPU 


def customiseTICLFromReco(process):
    # TensorFlow ESSource
    process.TFESSource = cms.Task(process.trackdnn_source)

    process.hgcalLayerClustersTask = cms.Task(process.hgcalLayerClustersEE,
                                              process.hgcalLayerClustersHSi,
                                              process.hgcalLayerClustersHSci,
                                              process.hgcalMergeLayerClusters)

# Reconstruction

    process.TICL = cms.Path(process.hgcalLayerClustersTask,
                            process.TFESSource,
                            process.ticlLayerTileTask,
                            process.ticlIterationsTask,
#                            process.ticlGraphTask,
                            process.ticlTracksterMergeTask,
                            process.ticlSimTrackstersTask)
# Validation
    process.TICL_ValidationProducers = cms.Task(process.hgcalRecHitMapProducer,
                                                process.lcAssocByEnergyScoreProducer,
                                                process.layerClusterCaloParticleAssociationProducer,
                                                process.scAssocByEnergyScoreProducer,
                                                process.layerClusterSimClusterAssociationProducer,
                                                process.simTsAssocByEnergyScoreProducer,
                                                process.simTracksterHitLCAssociatorByEnergyScoreProducer,
                                                process.tracksterSimTracksterAssociationLinking,
                                                process.tracksterSimTracksterAssociationPR,
                                                process.tracksterSimTracksterAssociationLinkingbyCLUE3D,
                                                process.tracksterSimTracksterAssociationPRbyCLUE3D,
                                                process.tracksterSimTracksterAssociationLinkingPU,
                                                process.tracksterSimTracksterAssociationPRPU,
  #                                              process.tracksterSimTracksterAssociationLinkingbyCLUE3DPU,
   #                                             process.tracksterSimTracksterAssociationPRbyCLUE3DPU 
                                                )
    process.TICL_Validator = cms.Task(process.hgcalValidator)
    process.TICL_Validation = cms.Path(process.TICL_ValidationProducers,
                                       process.TICL_Validator
                                       )
# Path and EndPath definitions
    process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
    process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
    process.schedule = cms.Schedule(process.TICL,
                                    process.TICL_Validation,
                                    process.FEVTDEBUGHLToutput_step,
                                    process.DQMoutput_step)
# call to customisation function customiseHGCalOnlyEventContent imported from RecoHGCal.Configuration.RecoHGCal_EventContent_cff
    process = customiseHGCalOnlyEventContent(process)

    return process


def customiseTICLForDumper(process):
    
    process.simpleValidation = simpleValidation.clone()

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
    )
    process.ticlGraphAnalyzer = ticlGraphAnalyzer.clone(
    )
    process.ticlGraphAnalyzerCone = ticlGraphAnalyzer.clone(
            ticlGraph = 'ticlGraph:cone' 
    )
    process.TFileService = cms.Service("TFileService",
                                       fileName=cms.string("histo.root")
                                       )
    process.FEVTDEBUGHLToutput_step = cms.EndPath(
    process.FEVTDEBUGHLToutput + process.simpleValidation + process.ticlDumper)
    #process.FEVTDEBUGHLToutput + process.ticlDumper + process.ticlGraphAnalyzer + process.ticlGraphAnalyzerCone)
    return process

def customiseTICLForValidationPlot(process):
    
    process.simpleValidation = simpleValidation.clone()
    process.FEVTDEBUGHLToutput = cms.EndPath(process.simpleValidation)

    return process
