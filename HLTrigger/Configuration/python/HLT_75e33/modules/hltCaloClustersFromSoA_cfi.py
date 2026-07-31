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

hltCaloClustersFromSoASerialSync = hltCaloClustersFromSoA.clone(
    src = "hltHgCalLayerClustersFromSoAProducerSerialSync"
)

from Configuration.ProcessModifiers.alpaka_cff import alpaka
from Configuration.ProcessModifiers.ticl_barrel_cff import ticl_barrel

(alpaka & ~ticl_barrel).toModify(hltCaloClustersFromSoA,
    src = "hltHgCalLayerClustersFromSoAProducer"
)
