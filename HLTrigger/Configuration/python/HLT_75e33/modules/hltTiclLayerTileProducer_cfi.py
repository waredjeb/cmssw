import FWCore.ParameterSet.Config as cms

hltTiclLayerTileProducer = cms.EDProducer("TICLLayerTileProducer",
    detector = cms.string('HGCAL'),
    layer_clusters = cms.InputTag("hltMergeLayerClusters"),
    mightGet = cms.optional.untracked.vstring
)

from Configuration.ProcessModifiers.alpaka_cff import alpaka
from Configuration.ProcessModifiers.ticl_barrel_cff import ticl_barrel

(alpaka & ~ticl_barrel).toModify(hltTiclLayerTileProducer,
    layer_clusters = "hltHgCalLayerClustersFromSoAProducer"
)
