import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.ticlLayerTileProducer_cfi import ticlLayerTileProducer
from RecoTracker.IterativeTracking.iterativeTk_cff import trackdnn_source
from RecoLocalCalo.HGCalRecProducers.hgcalRecHitMapProducer_cfi import hgcalRecHitMapProducer

from RecoHGCal.TICL.CLUE3DEM_cff import *
from RecoHGCal.TICL.CLUE3DHAD_cff import *
from RecoHGCal.TICL.SimTracksters_cff import *
from RecoHGCal.TICL.pfTICLProducer_cfi import pfTICLProducer as _pfTICLProducer

from RecoHGCal.TICL.ticlLayerTileProducer_cfi import ticlLayerTileProducer
from RecoHGCal.TICL.pfTICLProducer_cfi import pfTICLProducer as _pfTICLProducer
from RecoHGCal.TICL.tracksterSelectionTf_cfi import *

from RecoHGCal.TICL.tracksterLinksProducer_cfi import tracksterLinksProducer as _tracksterLinksProducer
from RecoHGCal.TICL.ticlCandidateProducer_cfi import ticlCandidateProducer as _ticlCandidateProducer
from RecoHGCal.Configuration.RecoHGCal_EventContent_cff import customiseForTICLv5EventContent
from RecoHGCal.TICL.iterativeTICL_cff import ticlIterLabels, ticlIterLabelsMerge
from RecoHGCal.TICL.ticlDumper_cfi import ticlDumper
from RecoHGCal.TICL.mergedTrackstersProducer_cfi import mergedTrackstersProducer as _mergedTrackstersProducer
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterAssociationLinkingbyCLUE3D as _tracksterSimTracksterAssociationLinkingbyCLUE3D
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterAssociationPRbyCLUE3D  as _tracksterSimTracksterAssociationPRbyCLUE3D 
from Validation.HGCalValidation.HGCalValidator_cff import hgcalValidatorv5 


def customiseForTICLv5(process, enableDumper = False):
    #v = [float(1.35227069312378), float(0.597984541122229), 0.357752062837009, 0.190213241768714, 0.8, 1.29023531258318, 3.99098405073771, 2.64272990702104, 0.873443419035198, 0.837962914766715, 2.96740409778711, 2.48365834717034, 4.09302356886537, 1.62460320930781, 0.290721459465582, 0.5, 0, 5, 1, 0.36248203427106, 0.898885015049942, 1.5, 0.649266352277089, 0.05, 0.005, 0.0184891635833641, 0.334653355634394, 0.127805653259902, 0.5, 1.06819759109575, 1.87334588417966, 2, 7, 6.46911383351854, 7, 0.001, 0.0001, 0.0001, 0.973040101647016, 4.10681076598546, 5, 5, 2.5287882522632, 4.05811730381901, 5, 0.05, 0.349565355817902, 0.374412328923012, 5, 3.52770641432807, 4.53713490741646, 0, 3, 3, 0.382320673078743, 0.5, 1.00419486137907, 1.8, 483.23572131575, 0.98, 10, 0.9, 400]
    v = [0.599999999999999978, 0.599999999999999978, 0.149999999999999994, 0.149999999999999994, 1.800000000000000044, 1.800000000000000044, 5.000000000000000000, 5.000000000000000000, 0.500000000000000000, 0.000000000000000000, 3.000000000000000000, 3.000000000000000000, 3.240000000000000213, 3.240000000000000213, 0.200000000000000011, 0.200000000000000011, 2.000000000000000000, 2.000000000000000000, 0.200000000000000011, 0.200000000000000011, 0.100000000000000006, 0.358835416872817725, 0.101798260544287922, 0.005000000000000000, 0.036037821509803454, 0.019492355642557312, 0.500000000000000000, 0.500000000000000000, 0.050000000000000003, 1.456226335177476461, 1.181350011595288407, 1.632590158037479533, 1.000000000000000000, 7.000000000000000000000000, 7.000000000000000000, 0.000858981662552112, 0.000902743678359972, 0.000100000000000000, 1.000000000000000000, 1.000000000000000000, 5.000000000000000000, 1.000000000000000000, 4.802121152124460224, 5.000000000000000000, 1.000000000000000000, 0.691276700585935155, 0.059240170020699423, 1.000000000000000000, 0.000000000000000000, 3.915556919195346097, 0.141372798471569361, 0.362057089490857908, 2.224109837037461901, 0.565955771581091671, 0.599999999999999978, 1.568465062992860481, 1.601329925978132618, 1.211976274452753355, 132.357016684851657828, 0.800000000000000044, 10.000000000000000000, 0.936285648615934107, 2954.376808821650683967, 0.998446822352248908, 0.000000000000000000, 0.315906557678854616, 0.251000000000000001, 0.249193075117370905, 0.000120197508546043]


    print(v[0])
    process.ticlLayerTileTask = cms.Task(ticlLayerTileProducer)
    process.ticlSeedingGlobal = ticlSeedingGlobal.copy()
    process.filteredLayerClustersCLUE3DEM = filteredLayerClustersCLUE3DEM.copy()
    process.ticlTrackstersCLUE3DEM = ticlTrackstersCLUE3DEM.copy()

    process.ticlCLUE3DEMStepTask = cms.Task(process.ticlSeedingGlobal
        ,process.filteredLayerClustersCLUE3DEM
        ,process.ticlTrackstersCLUE3DEM)

    process.filteredLayerClustersCLUE3DHAD = filteredLayerClustersCLUE3DHAD.copy()
    process.ticlTrackstersCLUE3DHAD = ticlTrackstersCLUE3DHAD.copy()

    process.ticlCLUE3DHADStepTask = cms.Task(process.ticlSeedingGlobal
        ,process.filteredLayerClustersCLUE3DHAD
        ,process.ticlTrackstersCLUE3DHAD)

    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.criticalDensity = [v[0], v[1], 0.6]
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.criticalEtaPhiDistance = [0.025, 0.025, 0.025]
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.criticalSelfDensity = [v[2], v[3], 0.15]
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.criticalXYDistance = [v[4],v[5], 1.8]
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.criticalZDistanceLyr = [int(v[6]) ,int(v[7]), 5]
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.cutHadProb = cms.double(v[8])
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.densityEtaPhiDistanceSqr = [0.0008, 0.0008, 0.0008]
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.densityOnSameLayer = cms.bool(bool(int(v[9])))
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.densitySiblingLayers = [int(v[10]), int(v[11]), 3]
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.densityXYDistanceSqr = [v[12],v[13],3.24]
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.doPidCut = cms.bool(True)
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.eid_input_name = cms.string('input')
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.eid_min_cluster_energy = cms.double(1)
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.eid_n_clusters = cms.int32(10)
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.eid_n_layers = cms.int32(50)
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.eid_output_name_energy = cms.string('output/regressed_energy')
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.eid_output_name_id = cms.string('output/id_probabilities')
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.kernelDensityFactor = [v[14], v[15], 0.2]
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.minNumLayerCluster = [int(v[16]),int(v[17]),2]
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.nearestHigherOnSameLayer = cms.bool(False)
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.outlierMultiplier = [v[18],v[19],0.2]
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.rescaleDensityByZ = cms.bool(False)
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.type = cms.string('CLUE3D')
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.useAbsoluteProjectiveScale = cms.bool(True)
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.useClusterDimensionXY = cms.bool(False)

    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.criticalDensity = cms.vdouble(v[20],v[21],v[22])
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.criticalEtaPhiDistance = cms.vdouble(v[23], v[24], v[25])
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.criticalSelfDensity = cms.vdouble(v[26],v[27],v[28])
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.criticalXYDistance = cms.vdouble(v[29], v[30], v[31])
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.criticalZDistanceLyr = cms.vint32(int(v[32]), int(v[33]), int(v[34]))
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.cutHadProb = cms.double(0.5)
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.densityEtaPhiDistanceSqr = cms.vdouble(v[35], v[36], v[37])
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.densityOnSameLayer = cms.bool(bool(int(v[38])))
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.densitySiblingLayers = cms.vint32(int(v[39]), int(v[40]), int(v[41]))
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.densityXYDistanceSqr = cms.vdouble(v[42], v[43], v[44])
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.doPidCut = cms.bool(False)
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.eid_input_name = cms.string('input')
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.eid_min_cluster_energy = cms.double(1)
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.eid_n_clusters = cms.int32(10)
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.eid_n_layers = cms.int32(50)
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.eid_output_name_energy = cms.string('output/regressed_energy')
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.eid_output_name_id = cms.string('output/id_probabilities')
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.kernelDensityFactor = cms.vdouble(v[45], v[46], v[47])
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.minNumLayerCluster = cms.vint32(int(v[48]), int(v[49]), int(v[50]))
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.nearestHigherOnSameLayer = cms.bool(False)
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.outlierMultiplier = cms.vdouble(v[51], v[52], v[53])
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.rescaleDensityByZ = cms.bool(False)
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.type = cms.string('CLUE3D')
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.useAbsoluteProjectiveScale = cms.bool(True)
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.useClusterDimensionXY = cms.bool(False)
    from SimTracker.TrackerHitAssociation.tpClusterProducerDefault_cfi import tpClusterProducerDefault as _tpClusterProducerDefault


    process.tpClusterProducer = _tpClusterProducerDefault.clone(
      trackingParticleSrc = "mixData:MergedTrackTruth",
      pixelSimLinkSrc = "mixData:PixelDigiSimLink",
      stripSimLinkSrc = "mixData:StripDigiSimLink",
      phase2OTSimLinkSrc = "mixData:Phase2OTDigiSimLink",
        
     )

    process.quickTrackAssociatorByHits = cms.EDProducer("QuickTrackAssociatorByHitsProducer",
      AbsoluteNumberOfHits = cms.bool(False),
      Cut_RecoToSim = cms.double(0.75),
      SimToRecoDenominator = cms.string('reco'), # either "sim" or "reco"
      Quality_SimToReco = cms.double(0.5),
      Purity_SimToReco = cms.double(0.75),
      ThreeHitTracksAreSpecial = cms.bool(True),
            PixelHitWeight = cms.double(1.0),
            useClusterTPAssociation = cms.bool(True),
            cluster2TPSrc = cms.InputTag("tpClusterProducer")
    )

    process.trackingParticleRecoTrackAsssociation = cms.EDProducer("TrackAssociatorEDProducer",
        associator = cms.InputTag('quickTrackAssociatorByHits'),
        label_tp = cms.InputTag("mix","MergedTrackTruth"),
        label_tr = cms.InputTag("generalTracks"),
        ignoremissingtrackcollection=cms.untracked.bool(False)
    )

    process.ticlIterationsTask = cms.Task(
        process.ticlCLUE3DEMStepTask,
        process.ticlCLUE3DHADStepTask,
    )

    #process.tpRecoAssociationTask = cms.Task(process.trackingParticleRecoTrackAsssociation)

    process.ticlTracksterLinks = _tracksterLinksProducer.clone()
    process.ticlTracksterLinks.linkingPSet.track_time_quality_threshold = cms.double(0.5)
    process.ticlTracksterLinks.linkingPSet.wind = cms.double(v[54])
    process.ticlTracksterLinks.linkingPSet.angle0 = cms.double(v[55])
    process.ticlTracksterLinks.linkingPSet.angle1 = cms.double(v[56])
    process.ticlTracksterLinks.linkingPSet.angle2 = cms.double(v[57])
    process.ticlTracksterLinks.linkingPSet.maxConeHeight = cms.double(v[58])
    process.ticlTracksterLinks.linkingPSet.pcaQuality = cms.double(v[59])
    process.ticlTracksterLinks.linkingPSet.pcaQualityLCSize = cms.uint32(int(v[60]))
    process.ticlTracksterLinks.linkingPSet.dotProdCut = cms.double(v[61])
    process.ticlTracksterLinks.linkingPSet.maxDistSkeletonsSq = cms.double(v[62])
    process.ticlTracksterLinks.linkingPSet.angle0_scaling = cms.double(v[63])
    process.ticlTracksterLinks.linkingPSet.angle1_scaling = cms.double(v[64])
    process.ticlTracksterLinks.linkingPSet.angle2_scaling = cms.double(v[65])
    process.ticlTracksterLinks.linkingPSet.algo_verbosity = cms.int32(0)
    process.ticlTracksterLinks.tracksters_collections = cms.VInputTag("ticlTrackstersCLUE3DEM", "ticlTrackstersCLUE3DHAD")

    process.ticlTracksterLinksTask = cms.Task(process.ticlTracksterLinks)

    process.ticlCandidate = _ticlCandidateProducer.clone()
    process.ticlCandidateTask = cms.Task(process.ticlCandidate)

    process.tracksterSimTracksterAssociationLinkingbyCLUE3DEM = _tracksterSimTracksterAssociationLinkingbyCLUE3D.clone(
        label_tst = cms.InputTag("ticlTrackstersCLUE3DEM")
        )
    process.tracksterSimTracksterAssociationPRbyCLUE3DEM = _tracksterSimTracksterAssociationPRbyCLUE3D.clone(
        label_tst = cms.InputTag("ticlTrackstersCLUE3DEM")
        )
    process.tracksterSimTracksterAssociationLinkingbyCLUE3DHAD = _tracksterSimTracksterAssociationLinkingbyCLUE3D.clone(
        label_tst = cms.InputTag("ticlTrackstersCLUE3DHAD")
        )
    process.tracksterSimTracksterAssociationPRbyCLUE3DHAD = _tracksterSimTracksterAssociationPRbyCLUE3D.clone(
        label_tst = cms.InputTag("ticlTrackstersCLUE3DHAD")
        )

    process.mergedTrackstersProducer = _mergedTrackstersProducer.clone()    

    process.tracksterSimTracksterAssociationLinkingbyCLUE3D = _tracksterSimTracksterAssociationLinkingbyCLUE3D.clone(
        label_tst = cms.InputTag("mergedTrackstersProducer")
        )
    process.tracksterSimTracksterAssociationPRbyCLUE3D = _tracksterSimTracksterAssociationPRbyCLUE3D.clone(
        label_tst = cms.InputTag("mergedTrackstersProducer")
        )

    process.TFESSource = cms.Task(process.trackdnn_source)
    process.iterTICLTask = cms.Task(process.ticlLayerTileTask,
                                     process.TFESSource,
                                     process.ticlIterationsTask,
                                     process.mergedTrackstersProducer,
                                     process.ticlTracksterLinksTask)
                                     #process.ticlCandidateTask)
    process.particleFlowClusterHGCal.initialClusteringStep.tracksterSrc = "ticlTracksterLinks"
    process.globalrecoTask.remove(process.ticlTrackstersMerge)

    process.tracksterSimTracksterAssociationLinking.label_tst = cms.InputTag("ticlTracksterLinks")
    process.tracksterSimTracksterAssociationPR.label_tst = cms.InputTag("ticlTracksterLinks")

    process.tracksterSimTracksterAssociationLinkingPU.label_tst = cms.InputTag("ticlTracksterLinks")
    process.tracksterSimTracksterAssociationPRPU.label_tst = cms.InputTag("ticlTracksterLinks")
    process.mergeTICLTask = cms.Task()
    #process.pfTICL.ticlCandidateSrc = cms.InputTag("ticlCandidate") 
    process.hgcalAssociators = cms.Task(process.mergedTrackstersProducer,process.hgcalRecHitMapProducer,
                            process.layerClusterCaloParticleAssociationProducer,
                            process.lcAssocByEnergyScoreProducer,
                            process.scAssocByEnergyScoreProducer, process.layerClusterSimClusterAssociationProducer,
                            process.lcSimTSAssocByEnergyScoreProducer, process.layerClusterSimTracksterAssociationProducer,
                            process.simTsAssocByEnergyScoreProducer,  process.simTracksterHitLCAssociatorByEnergyScoreProducer,
                            process.tracksterSimTracksterAssociationLinking, process.tracksterSimTracksterAssociationPR,
                            process.tracksterSimTracksterAssociationLinkingbyCLUE3D, process.tracksterSimTracksterAssociationPRbyCLUE3D,
                            process.tracksterSimTracksterAssociationLinkingbyCLUE3DEM, process.tracksterSimTracksterAssociationPRbyCLUE3DEM,
                            process.tracksterSimTracksterAssociationLinkingbyCLUE3DHAD, process.tracksterSimTracksterAssociationPRbyCLUE3DHAD,
                            process.tracksterSimTracksterAssociationLinkingPU, process.tracksterSimTracksterAssociationPRPU
                            )

    process.globalPrevalidationHGCal = cms.Sequence(process.hgcalAssociators) 
    process.hgcalAssociatorsPath = cms.Path(process.globalPrevalidationHGCal)
    
    process.hgcalValidatorv5 = hgcalValidatorv5.clone(
        ticlTrackstersMerge = cms.InputTag("ticlTracksterLinks"),
        trackstersclue3d = cms.InputTag("mergedTrackstersProducer")
    )
    process.hgcalValidatorSequence = cms.Sequence(process.hgcalValidatorv5)
    process.hgcalValidation = cms.Sequence(process.hgcalSimHitValidationEE+process.hgcalSimHitValidationHEF+process.hgcalSimHitValidationHEB+process.hgcalDigiValidationEE+process.hgcalDigiValidationHEF+process.hgcalDigiValidationHEB+process.hgcalRecHitValidationEE+process.hgcalRecHitValidationHEF+process.hgcalRecHitValidationHEB+process.hgcalHitValidationSequence+process.hgcalValidatorSequence+process.hgcalTiclPFValidation+process.hgcalPFJetValidation)
    process.globalValidationHGCal = cms.Sequence(process.hgcalValidation)
    process.validation_step9 = cms.EndPath(process.globalValidationHGCal)
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
            trackstersclue3d = cms.InputTag('mergedTrackstersProducer'),
            #ticlcandidates = cms.InputTag("ticlCandidate"),
            trackstersmerged = cms.InputTag("ticlTracksterLinks")
        )
        process.TFileService = cms.Service("TFileService",
                                           fileName=cms.string("histoBestConfigTICLv5.root")
                                           )

        process.TICL = cms.Path(process.iterTICLTask)
        process.FEVTDEBUGHLToutput_step = cms.EndPath(
            process.FEVTDEBUGHLToutput + process.ticlDumper)



        # Schedule definition
        process.schedule = cms.Schedule(process.TICL,
                                        process.hgcalAssociatorsPath,
                                        process.FEVTDEBUGHLToutput_step)
            #                            process.DQMoutput_step)
    process = customiseForTICLv5EventContent(process)

    return process
