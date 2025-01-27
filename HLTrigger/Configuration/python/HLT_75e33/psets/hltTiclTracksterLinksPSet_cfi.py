import FWCore.ParameterSet.Config as cms

hltTiclTracksterLinksBySkeletonsPSet = cms.PSet(
      cylinder_radius_sqr_split = cms.double(9),
      proj_distance_split = cms.double(5),
      track_time_quality_threshold = cms.double(0.5),
      min_num_lcs = cms.uint32(15),
      min_trackster_energy = cms.double(20),
      pca_quality_th = cms.double(0.85),
      dot_prod_th = cms.double(0.97),
      deltaRxy = cms.double(4),
      lower_boundary = cms.vdouble(
        20,
        10
      ),
      upper_boundary = cms.vdouble(
        150,
        100
      ),
      upper_distance_projective_sqr = cms.vdouble(
        30,
        60
      ),
      lower_distance_projective_sqr = cms.vdouble(
        30,
        60
      ),
      min_distance_z = cms.vdouble(
        35,
        35
      ),
      upper_distance_projective_sqr_closest_points = cms.vdouble(
        5,
        30
      ),
      lower_distance_projective_sqr_closest_points = cms.vdouble(
        10,
        50
      ),
      max_z_distance_closest_points = cms.vdouble(
        35,
        35
      ),
      cylinder_radius_sqr = cms.vdouble(
        9,
        15
      ),
      algo_verbosity = cms.int32(0),
      type = cms.string('Skeletons')
    )

hltTiclTracksterLinksBySuperClusteringDNNPSet = cms.PSet(
      algo_verbosity = cms.int32(0),
      onnxModelPath = cms.FileInPath('RecoHGCal/TICL/data/superclustering/supercls_v2p1.onnx'),
      dnnInputsVersion = cms.string('v2'),
      inferenceBatchSize = cms.uint32(100000),
      nnWorkingPoint = cms.double(0.3),
      deltaEtaWindow = cms.double(0.1),
      deltaPhiWindow = cms.double(0.5),
      seedPtThreshold = cms.double(4),
      candidateEnergyThreshold = cms.double(2),
      explVarRatioCut_energyBoundary = cms.double(50),
      explVarRatioMinimum_lowEnergy = cms.double(0.92),
      explVarRatioMinimum_highEnergy = cms.double(0.95),
      filterByTracksterPID = cms.bool(True),
      tracksterPIDCategoriesToFilter = cms.vint32(
        0,
        1
      ),
      PIDThreshold = cms.double(0.8),
      type = cms.string('SuperClusteringDNN')
    )

hltTiclTracksterLinksBySuperClusteringMustachePSet = cms.PSet(
      algo_verbosity = cms.int32(0),
      seedThresholdPt = cms.double(1),
      candidateEnergyThreshold = cms.double(0.15),
      filterByTracksterPID = cms.bool(True),
      tracksterPIDCategoriesToFilter = cms.vint32(
        0,
        1
      ),
      PIDThreshold = cms.double(0.8),
      type = cms.string('SuperClusteringMustache')
    )
