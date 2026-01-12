import FWCore.ParameterSet.Config as cms
from RecoHGCal.TICL.tracksterCleaningProducer_cfi import tracksterCleaningProducer as _tracksterCleaningProducer

hltTiclTracksterCleaning = _tracksterCleaningProducer.clone(
    linkedTracksters      = cms.InputTag("hltTiclTracksterLinks"),
    clue3DTracksters      = cms.InputTag("hltTiclTrackstersCLUE3DHigh"),
    clue3DInLinkedIndices = cms.InputTag("hltTiclTracksterLinks","linkedTracksterIdToInputTracksterId"),
    algo_verbosity = cms.int32(0),
    cleaner = cms.PSet(
      type = cms.string('Beta'),
      algo_verbosity = cms.int32(0),
      betaContamMin = cms.double(1.5),
      R0 = cms.double(0.1),
      useRawEnergy = cms.bool(True),
      epsE = cms.double(1e-06),
      epsDR = cms.double(1e-06),
      weightMode = cms.bool(True),
      emitDroppedAsStandalone = cms.bool(False),
      zAbsCut = cms.double(25),
      tAbsCut = cms.double(0.15),
      sigmaZ = cms.double(12.5),
      sigmaT = cms.double(0.08),
      sigmaDR = cms.double(0.08),
      zPower = cms.double(1.5),
      tPower = cms.double(0.5),
      drPower = cms.double(0.5),
      wmin = cms.double(0.001),
      doPruning = cms.bool(False),
      pruneWmin = cms.double(0.01),
      pruneUseSeparateKernels = cms.bool(False),
      sigmaZ_prune = cms.double(12.5),
      sigmaT_prune = cms.double(0.08),
      sigmaDR_prune = cms.double(0.08),
      zPower_prune = cms.double(1.5),
      tPower_prune = cms.double(0.5),
      drPower_prune = cms.double(0.5)
    ),
)
