import FWCore.ParameterSet.Config as cms

hltFilteredLayerClustersRecovery = cms.EDProducer("FilteredLayerClustersProducer",
    LayerClusters = cms.InputTag("hltMergeLayerClusters"),
    LayerClustersInputMask = cms.InputTag("hltTiclTrackstersCLUE3DHigh"),
    algo_number = cms.vint32(6, 7, 8),
    clusterFilter = cms.string('ClusterFilterByAlgoAndSize'),
    iteration_label = cms.string('Recovery'),
    max_cluster_size = cms.int32(9999),
    max_layerId = cms.int32(9999),
    mightGet = cms.optional.untracked.vstring,
    min_cluster_size = cms.int32(2),
    min_layerId = cms.int32(0)
)

from Configuration.ProcessModifiers.alpaka_cff import alpaka
from Configuration.ProcessModifiers.ticl_barrel_cff import ticl_barrel

(alpaka & ~ticl_barrel).toModify(hltFilteredLayerClustersRecovery,
    LayerClusters = "hltHgCalLayerClustersFromSoAProducer"
)
