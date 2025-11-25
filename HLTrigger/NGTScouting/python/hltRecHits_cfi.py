import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.nano_cff import nanoMetadata
from Validation.HGCalValidation.HLT_TICLIterLabels_cff import hltTiclIterLabels

# Tracksters
hltRecHitsTable = cms.EDProducer("RecHitsExtraTableProducer",
    tableName=cms.string("HGCalRecHits"),
    skipNonExistingSrc=cms.bool(False),
    recHits=cms.VInputTag([("hltHGCalRecHit", "HGCEERecHits"),("hltHGCalRecHit", "HGCHEBRecHits"), ("hltHGCalRecHit", "HGCHEFRecHits")]),
    precision=cms.int32(7))


hltRecHitsTableSequence = cms.Sequence(hltRecHitsTable) 

