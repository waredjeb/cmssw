import FWCore.ParameterSet.Config as cms

# Barrel layer clustering (CLUEstering) for EB and HB.
#
# Each detector runs the same chain:
#   hltBarrelRecHitsSoA<Det>        (legacy CaloRecHit -> CaloRecHit SoA)
#   hltBarrelPFRecHitsSoA<Det>      (CaloRecHit SoA -> PFRecHit SoA)
#   hltBarrelCLUEstering<Det>       (CLUEstering clustering)
#   hltBarrelSoALayerClusters<Det>  (build the CaloCluster SoA)
#   hltBarrelLayerClusters<Det>     (SoA -> legacy reco::CaloCluster)
#
from ..modules.hltBarrelRecHitsSoAEB_cfi import *
from ..modules.hltBarrelPFRecHitsSoAEB_cfi import *
from ..modules.hltBarrelCLUEsteringEB_cfi import *
from ..modules.hltBarrelSoALayerClustersEB_cfi import *
from ..modules.hltBarrelLayerClustersEB_cfi import *
from ..modules.hltBarrelRecHitsSoAHB_cfi import *
from ..modules.hltBarrelPFRecHitsSoAHB_cfi import *
from ..modules.hltBarrelCLUEsteringHB_cfi import *
from ..modules.hltBarrelSoALayerClustersHB_cfi import *
from ..modules.hltBarrelLayerClustersHB_cfi import *

HLTBarrelLayerClustersSequence = cms.Sequence(
        hltBarrelRecHitsSoAEB+
        hltBarrelPFRecHitsSoAEB+
        hltBarrelCLUEsteringEB+
        hltBarrelSoALayerClustersEB+
        hltBarrelLayerClustersEB+
        hltBarrelRecHitsSoAHB+
        hltBarrelPFRecHitsSoAHB+
        hltBarrelCLUEsteringHB+
        hltBarrelSoALayerClustersHB+
        hltBarrelLayerClustersHB)
