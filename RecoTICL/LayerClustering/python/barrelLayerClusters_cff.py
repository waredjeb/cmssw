import FWCore.ParameterSet.Config as cms

#####################################################################
# CLUEstering barrel layer-clustering chain
#
# The chain per detector is:
#   barrelRecHitsSoA<Det>        (legacy CaloRecHit -> CaloRecHit SoA)
#   barrelPFRecHitsSoA<Det>      (CaloRecHit SoA -> PFRecHit SoA)
#   barrelCLUEstering<Det>       (CLUE clustering on the PFRecHit SoA)
#   barrelSoALayerClusters<Det>  (build the CaloCluster SoA)
#   barrelLayerClusters<Det>     (SoA -> legacy reco::CaloCluster products)
#
# The HCAL topology and parameters come from the 'offline' PF ESProducers
# (alpaka_pfClusteringHBHEHF_esProducersTask in particleFlowCluster_cff).
#####################################################################
from RecoParticleFlow.PFRecHitProducer.ecalRecHitSoAProducer_cfi import ecalRecHitSoAProducer
from RecoParticleFlow.PFRecHitProducer.hcalRecHitSoAProducer_cfi import hcalRecHitSoAProducer
from RecoParticleFlow.PFRecHitProducer.pfRecHitECALParamsESProducer_cfi import pfRecHitECALParamsESProducer
from RecoParticleFlow.PFRecHitProducer.pfRecHitECALTopologyESProducer_cfi import pfRecHitECALTopologyESProducer
from RecoParticleFlow.PFRecHitProducer.pfRecHitSoAProducerECAL_cfi import pfRecHitSoAProducerECAL
from RecoParticleFlow.PFRecHitProducer.pfRecHitSoAProducerHCAL_cfi import pfRecHitSoAProducerHCAL
from RecoTICL.LayerClustering.barrelCLUEsteringLayerClustersProducerECAL_cfi import barrelCLUEsteringLayerClustersProducerECAL
from RecoTICL.LayerClustering.barrelCLUEsteringLayerClustersProducerHCAL_cfi import barrelCLUEsteringLayerClustersProducerHCAL
from RecoTICL.LayerClustering.barrelSoALayerClustersProducerECAL_cfi import barrelSoALayerClustersProducerECAL
from RecoTICL.LayerClustering.barrelSoALayerClustersProducerHCAL_cfi import barrelSoALayerClustersProducerHCAL
from RecoTICL.LayerClustering.barrelLayerClustersFromSoAProducer_cfi import barrelLayerClustersFromSoAProducer

# ---- EB (ECAL barrel) ----
barrelPFRecHitECALTopologyRecordSource = cms.ESSource('EmptyESSource',
    recordName = cms.string('PFRecHitECALTopologyRecord'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)
# The PF ECAL thresholds are used both as PFRecHit energy cuts and as CLUE noise
barrelPFRecHitECALParamsESProducer = pfRecHitECALParamsESProducer.clone(appendToDataLabel = 'barrel')
barrelPFRecHitECALTopologyESProducer = pfRecHitECALTopologyESProducer.clone(appendToDataLabel = 'barrel')

barrelRecHitsSoAEB = ecalRecHitSoAProducer.clone(src = 'ecalRecHit:EcalRecHitsEB')
barrelPFRecHitsSoAEB = pfRecHitSoAProducerECAL.clone(
    producers = [cms.PSet(
        src = cms.InputTag('barrelRecHitsSoAEB'),
        params = cms.ESInputTag('barrelPFRecHitECALParamsESProducer', 'barrel'),
    )],
    topology = cms.ESInputTag('barrelPFRecHitECALTopologyESProducer', 'barrel'))
barrelCLUEsteringEB = barrelCLUEsteringLayerClustersProducerECAL.clone(
    recHits = 'barrelPFRecHitsSoAEB',
    noise = cms.ESInputTag('barrelPFRecHitECALParamsESProducer', 'barrel'),
    deltac = 1.8 * 0.0175,
    kappa = 1.)
barrelSoALayerClustersEB = barrelSoALayerClustersProducerECAL.clone(
    recHits = 'barrelPFRecHitsSoAEB',
    barrelRecHitsLayerClustersSoA = 'barrelCLUEsteringEB')
barrelLayerClustersEB = barrelLayerClustersFromSoAProducer.clone(
    src = 'barrelSoALayerClustersEB',
    recHits = 'barrelPFRecHitsSoAEB',
    barrelRecHitsLayerClustersSoA = 'barrelCLUEsteringEB',
    detector = 'EB')

# ---- HB (HCAL barrel) ----
barrelRecHitsSoAHB = hcalRecHitSoAProducer.clone(src = 'hbhereco')
barrelPFRecHitsSoAHB = pfRecHitSoAProducerHCAL.clone(
    producers = [cms.PSet(
        src = cms.InputTag('barrelRecHitsSoAHB'),
        params = cms.ESInputTag('pfRecHitHCALParamsESProducer', 'offline'),
    )],
    topology = cms.ESInputTag('pfRecHitHCALTopologyESProducer', 'offline'))
barrelCLUEsteringHB = barrelCLUEsteringLayerClustersProducerHCAL.clone(
    recHits = 'barrelPFRecHitsSoAHB',
    noise = cms.ESInputTag('pfRecHitHCALTopologyESProducer', 'offline'),
    deltac = 3 * 0.087,
    kappa = 0.)
barrelSoALayerClustersHB = barrelSoALayerClustersProducerHCAL.clone(
    recHits = 'barrelPFRecHitsSoAHB',
    barrelRecHitsLayerClustersSoA = 'barrelCLUEsteringHB')
barrelLayerClustersHB = barrelLayerClustersFromSoAProducer.clone(
    src = 'barrelSoALayerClustersHB',
    recHits = 'barrelPFRecHitsSoAHB',
    barrelRecHitsLayerClustersSoA = 'barrelCLUEsteringHB',
    detector = 'HB')

barrelLayerClustersEBTask = cms.Task(
    barrelPFRecHitECALTopologyRecordSource,
    barrelPFRecHitECALParamsESProducer,
    barrelPFRecHitECALTopologyESProducer,
    barrelRecHitsSoAEB, barrelPFRecHitsSoAEB, barrelCLUEsteringEB, barrelSoALayerClustersEB,
    barrelLayerClustersEB,
)
barrelLayerClustersHBTask = cms.Task(
    barrelRecHitsSoAHB, barrelPFRecHitsSoAHB, barrelCLUEsteringHB, barrelSoALayerClustersHB,
    barrelLayerClustersHB,
)
