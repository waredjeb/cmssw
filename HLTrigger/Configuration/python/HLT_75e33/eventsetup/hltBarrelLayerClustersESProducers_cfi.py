import FWCore.ParameterSet.Config as cms
from Configuration.ProcessModifiers.ticl_barrel_cff import ticl_barrel

# PF topology and thresholds for the alpaka barrel layer clustering; attached
# only under the ticl_barrel process modifier, since only that path consumes them
def _addHltBarrelLayerClustersESProducers(process):
    process.hltESSPFRecHitECALTopologyRecord = cms.ESSource('EmptyESSource',
        recordName = cms.string('PFRecHitECALTopologyRecord'),
        iovIsRunNotTime = cms.bool(True),
        firstValid = cms.vuint32(1)
    )
    process.hltESSPFRecHitHCALParamsRecord = cms.ESSource('EmptyESSource',
        recordName = cms.string('PFRecHitHCALParamsRecord'),
        iovIsRunNotTime = cms.bool(True),
        firstValid = cms.vuint32(1)
    )
    process.hltESSPFRecHitHCALTopologyRecord = cms.ESSource('EmptyESSource',
        recordName = cms.string('PFRecHitHCALTopologyRecord'),
        iovIsRunNotTime = cms.bool(True),
        firstValid = cms.vuint32(1)
    )
    # The PF ECAL thresholds are used both as PFRecHit energy cuts and as CLUE noise
    process.hltESPPFRecHitECALParams = cms.ESProducer('PFRecHitECALParamsESProducer@alpaka',
        cleaningThreshold = cms.double(2.0),
        appendToDataLabel = cms.string('')
    )
    process.hltESPPFRecHitECALTopology = cms.ESProducer('PFRecHitECALTopologyESProducer@alpaka',
        usePFThresholdsFromDB = cms.bool(False),
        appendToDataLabel = cms.string('')
    )
    process.hltESPPFRecHitHCALParams = cms.ESProducer('PFRecHitHCALParamsESProducer@alpaka',
        energyThresholdsHB = cms.vdouble(0.1, 0.2, 0.3, 0.3),
        energyThresholdsHE = cms.vdouble(0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
        appendToDataLabel = cms.string('')
    )
    # usePFThresholdsFromDB fills the seed thresholds used as CLUE noise
    process.hltESPPFRecHitHCALTopology = cms.ESProducer('PFRecHitHCALTopologyESProducer@alpaka',
        usePFThresholdsFromDB = cms.bool(True),
        appendToDataLabel = cms.string('')
    )

addHltBarrelLayerClustersESProducers = ticl_barrel.makeProcessModifier(_addHltBarrelLayerClustersESProducers)
