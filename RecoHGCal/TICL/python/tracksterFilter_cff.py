import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.tracksterFilterProducer_cfi import tracksterFilterProducer as _tracksterFilterProducer

# Filter to keep only EM tracksters (filter out hadronic)
filteredTrackstersEM = _tracksterFilterProducer.clone(
    tracksters = cms.InputTag("ticlTrackstersCLUE3DHigh"),
    trackstersMask = cms.InputTag("ticlTrackstersCLUE3DHigh", "tracksterMask"),
    filterParams = cms.PSet(
        keepHadronic = cms.bool(False)
    )
)

tracksterFilterEMTask = cms.Task(filteredTrackstersEM)
