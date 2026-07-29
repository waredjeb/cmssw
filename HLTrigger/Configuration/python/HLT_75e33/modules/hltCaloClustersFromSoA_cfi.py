import FWCore.ParameterSet.Config as cms

# Legacy AoS view of the merged layer clusters, for the edm::Ref consumers
# (associators, validation, PF, EGamma) that cannot read the SoA.
hltCaloClustersFromSoA = cms.EDProducer("CaloClustersFromSoAProducer",
    src = cms.InputTag("hltMergeLayerClusters"),
    timeClname = cms.string("timeLayerCluster"),
    mightGet = cms.optional.untracked.vstring
)

hltCaloClustersFromSoAL1Seeded = hltCaloClustersFromSoA.clone(
    src = "hltMergeLayerClustersL1Seeded"
)

# Used by the alpakaValidationHLT GPU-vs-CPU comparison.
hltCaloClustersFromSoASerialSync = hltCaloClustersFromSoA.clone(
    src = "hltMergeLayerClustersSerialSync"
)
