import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var,CandVars
from DPGAnalysis.HGCalNanoAOD.trackstersNanoAOD_cfi import SimTrackstersTables, TrackstersTable  

nanoMetadata = cms.EDProducer("UniqueStringProducer",
    strings = cms.PSet(
        tag = cms.string("RecoHGCal"),
    )
)

namedProducers = []
for i, _producer in enumerate(TrackstersTable):
    label = f"ticlTracksterTable{i}"
    globals()[label] = _producer.clone()
    namedProducers.append(globals()[label])

# sim tracksters
for j, _producer in enumerate(SimTrackstersTables):
    label = f"ticlSimTracksterTable{j}"
    globals()[label] = _producer.clone()
    namedProducers.append(globals()[label])

trackstersSeq = cms.Sequence(sum(namedProducers, cms.Sequence()))
trackstersNanoAODTask = cms.Task(nanoMetadata, *namedProducers)

