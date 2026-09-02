import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal, ticlSeedingGlobalHFNose
from RecoHGCal.TICL.trackstersProducer_cfi import trackstersProducer as _trackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer

# CLUSTER FILTERING/MASKING

filteredLayerClustersCLUE3DHigh = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSize",
    min_cluster_size = 2, # inclusive
    iteration_label = "CLUE3DHigh"
)

# PATTERN RECOGNITION

ticlTrackstersCLUE3DHigh = _trackstersProducer.clone(
    filtered_mask = "filteredLayerClustersCLUE3DHigh:CLUE3DHigh",
    seeding_regions = "ticlSeedingGlobal",
    itername = "CLUE3DHigh",
    patternRecognitionBy = "CLUE3D",
    pluginPatternRecognitionByCLUE3D = dict (
        criticalDensity = [0.6, 0.6, 0.6],
        criticalEtaPhiDistance = [0.025, 0.025, 0.025],
        kernelDensityFactor = [0.2, 0.2, 0.2],
        algo_verbosity = 0,
        doPidCut = True,
        cutHadProb = 999
    ),
    inferenceAlgo = cms.string('TracksterInferenceByCNN'),
    pluginInferenceAlgoTracksterInferenceByCNN = cms.PSet(
        algo_verbosity = cms.int32(0),
        type = cms.string("TracksterInferenceByCNN"),
        onnxModelPath = cms.string("RecoHGCal/TICL/data/ticlv5/onnx_models/CNN/patternrecognition/id_v0.onnx"),
        inputNames = cms.vstring("input"),
        outputNames = cms.vstring("pid_output"),
        eid_min_cluster_energy = cms.double(1.0),
        eid_n_layers = cms.int32(50),
        eid_n_clusters = cms.int32(10),
        doPID = cms.int32(1),
        miniBatchSize = cms.untracked.int32(64),
    ),
    pluginInferenceAlgoTracksterInferenceByDNN = cms.PSet(
        algo_verbosity = cms.int32(0),
        onnxPIDModelPath = cms.string('RecoHGCal/TICL/data/ticlv5/onnx_models/DNN/patternrecognition/id_v0.onnx'),
        onnxEnergyModelPath = cms.string(''),
        inputNames  = cms.vstring('input'),
        output_en   = cms.vstring('enreg_output'),
        output_id   = cms.vstring('pid_output'),
        eid_min_cluster_energy = cms.double(1),
        eid_n_layers = cms.int32(50),
        eid_n_clusters = cms.int32(10),
        doPID = cms.int32(1),
        doRegression = cms.int32(0),
        type = cms.string('TracksterInferenceByDNN')
    ),

    pluginInferenceAlgoTracksterInferenceByPFN = cms.PSet(
        algo_verbosity = cms.int32(0),
        onnxPIDModelPath = cms.string('RecoHGCal/TICL/data/ticlv5/onnx_models/PFN/patternrecognition/id_v0.onnx'),
        onnxEnergyModelPath = cms.string(''),
        inputNames  = cms.vstring('input','input_tr_features'),
        output_en   = cms.vstring('enreg_output'),
        output_id   = cms.vstring('pid_output'),
        eid_min_cluster_energy = cms.double(1),
        eid_n_layers = cms.int32(50),
        eid_n_clusters = cms.int32(10),
        doPID = cms.int32(1),
        doRegression = cms.int32(0),
        type = cms.string('TracksterInferenceByPFN')
    )


)

ticlCLUE3DHighStepTask = cms.Task(ticlSeedingGlobal
    ,filteredLayerClustersCLUE3DHigh
    ,ticlTrackstersCLUE3DHigh)


# with ticl_dev 3D pattern recognition swapped to CLUEstering.
from Configuration.ProcessModifiers.ticl_dev import ticl_dev
from RecoHGCal.TICL.trackstersCLUEsteringProducer_cfi import trackstersCLUEsteringProducer as _trackstersCLUEsteringProducer

ticlTrackstersCLUEsteringAssignment = _trackstersCLUEsteringProducer.clone()

ticl_dev.toModify(ticlTrackstersCLUE3DHigh,
    patternRecognitionBy = 'CLUEstering',
    pluginPatternRecognitionByCLUEstering = dict(
        tracksterAssignment = 'ticlTrackstersCLUEsteringAssignment',
        algo_verbosity = 0,
        doPidCut = True,
        cutHadProb = 999,
    )
)

ticl_dev.toReplaceWith(ticlCLUE3DHighStepTask, cms.Task(ticlSeedingGlobal
    ,filteredLayerClustersCLUE3DHigh
    ,ticlTrackstersCLUEsteringAssignment
    ,ticlTrackstersCLUE3DHigh))

