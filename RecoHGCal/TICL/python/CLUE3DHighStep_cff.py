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
    pluginInferenceAlgoTracksterInferenceByDNN = cms.PSet(
      algo_verbosity = cms.int32(0),
      onnxPIDModelPath = cms.FileInPath('RecoHGCal/TICL/data/RecoHGCal-TICL/ticlv5/onnx_models/patternrecognition/id_v0.onnx'),
      onnxEnergyModelPath = cms.FileInPath('RecoHGCal/TICL/data/RecoHGCal-TICL/ticlv5/onnx_models/patternrecognition/energy_v0.onnx'),
      eid_min_cluster_energy = cms.double(1),
      eid_n_layers = cms.int32(50),
      eid_n_clusters = cms.int32(10),
      doPID = cms.int32(1),
      doRegression = cms.int32(0),
      type = cms.string('TracksterInferenceByDNN')
    ),
    pluginInferenceAlgoTracksterInferenceByANN = cms.PSet(
      algo_verbosity = cms.int32(0),
      type = cms.string('TracksterInferenceByANN')
    
    ),


)

from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
ticl_v5.toModify(ticlTrackstersCLUE3DHigh.pluginPatternRecognitionByCLUE3D, computeLocalTime = cms.bool(True))

ticlCLUE3DHighStepTask = cms.Task(ticlSeedingGlobal
    ,filteredLayerClustersCLUE3DHigh
    ,ticlTrackstersCLUE3DHigh)

