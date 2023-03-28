import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.FastJetStep_cff import *
from RecoHGCal.TICL.CLUE3DHighStep_cff import *
from RecoHGCal.TICL.CLUE3DLowStep_cff import *
from RecoHGCal.TICL.MIPStep_cff import *
from RecoHGCal.TICL.TrkEMStep_cff import *
from RecoHGCal.TICL.TrkStep_cff import *
from RecoHGCal.TICL.EMStep_cff import *
from RecoHGCal.TICL.HADStep_cff import *

from RecoHGCal.TICL.ticlLayerTileProducer_cfi import ticlLayerTileProducer
from RecoHGCal.TICL.pfTICLProducer_cfi import pfTICLProducer as _pfTICLProducer
from RecoHGCal.TICL.trackstersMergeProducer_cfi import trackstersMergeProducer as _trackstersMergeProducer
from RecoHGCal.TICL.trackstersMergeProducerV3_cfi import trackstersMergeProducerV3 as _trackstersMergeProducerV3
from RecoHGCal.TICL.ticlGraphProducer_cfi import ticlGraphProducer as _ticlGraphProducer
from RecoHGCal.TICL.tracksterSelectionTf_cfi import *

ticlLayerTileTask = cms.Task(ticlLayerTileProducer)

ticlTrackstersMerge = _trackstersMergeProducer.clone(
  linkingPSet = cms.PSet(
    cutTk = cms.string('1.48 < abs(eta) < 3.0 && pt > 1. && quality("highPurity") && hitPattern().numberOfLostHits("MISSING_OUTER_HITS") < 5'),
    delta_tk_ts_layer1 = cms.double(0.02),
    delta_tk_ts_interface = cms.double(0.03),
    delta_ts_em_had = cms.double(0.03),
    delta_ts_had_had = cms.double(0.03),
    track_time_quality_threshold = cms.double(0.5),
    algo_verbosity = cms.int32(0),
    type = cms.string('LinkingAlgoByDirectionGeometric')
  
  )
)
ticlTrackstersMergeV3 = _trackstersMergeProducerV3.clone()
ticlGraph = _ticlGraphProducer.clone()


pfTICL = _pfTICLProducer.clone()
ticlPFTask = cms.Task(pfTICL)

ticlIterationsTask = cms.Task(
    ticlCLUE3DHighStepTask
)

from Configuration.ProcessModifiers.clue3D_cff import clue3D
clue3D.toModify(ticlIterationsTask, func=lambda x : x.add(ticlCLUE3DHighStepTask,ticlCLUE3DLowStepTask))

from Configuration.ProcessModifiers.fastJetTICL_cff import fastJetTICL
fastJetTICL.toModify(ticlIterationsTask, func=lambda x : x.add(ticlFastJetStepTask))

from Configuration.ProcessModifiers.ticl_v3_cff import ticl_v3
ticl_v3.toModify(ticlIterationsTask, func=lambda x : x.add( ticlTrkEMStepTask
    ,ticlEMStepTask
    ,ticlTrkStepTask
    ,ticlHADStepTask) )
ticlIterLabels = [_step.itername.value() for _iteration in ticlIterationsTask for _step in _iteration if (_step._TypedParameterizable__type == "TrackstersProducer")]

ticlTracksterMergeTask = cms.Task(ticlTrackstersMerge)
ticlTracksterMergeTaskV3 = cms.Task(ticlTrackstersMergeV3)
ticlGraphTask = cms.Task(ticlGraph)

ticl_v3.toModify(pfTICL, ticlCandidateSrc = "ticlTrackstersMergeV3")

mergeTICLTask = cms.Task(ticlLayerTileTask
    ,ticlIterationsTask
    ,ticlTracksterMergeTask
    ,ticlPFTask
    ,ticlGraphTask
)

ticl_v3.toModify(mergeTICLTask, func=lambda x : x.add(ticlTracksterMergeTaskV3))
ticlIterLabelsMerge = ticlIterLabels + ["Merge"]

ticlIterLabelsMergeV3 = ticlIterLabels + ["MergeV3"]
ticl_v3.toModify(ticlIterLabelsMerge, func=lambda x : x.extend(ticlIterLabelsMergeV3))

iterTICLTask = cms.Task(mergeTICLTask
    ,ticlPFTask)

ticlLayerTileHFNose = ticlLayerTileProducer.clone(
    detector = 'HFNose'
)

ticlLayerTileHFNoseTask = cms.Task(ticlLayerTileHFNose)

iterHFNoseTICLTask = cms.Task(ticlLayerTileHFNoseTask
    ,ticlHFNoseTrkEMStepTask
    ,ticlHFNoseEMStepTask
    ,ticlHFNoseTrkStepTask
    ,ticlHFNoseHADStepTask
    ,ticlHFNoseMIPStepTask
)
