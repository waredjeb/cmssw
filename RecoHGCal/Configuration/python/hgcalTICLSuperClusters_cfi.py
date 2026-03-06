import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.nano_cff import nanoMetadata
from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5

hgcalUpgradeNanoTask = cms.Task(nanoMetadata)

# Offline TICL SuperClusters table - adapted from HLT version
# NOTE: The actual superclustering collection name in offline RECO is typically:
# - particleFlowSuperClusterHGCal (for PF-based mustache superclustering)
# - or output from ticlEGammaSuperClusterProducer (for TICL-based superclustering)
# Users should verify which superclustering is run in their workflow
ticlSuperClustersTable = cms.EDProducer(
    "TICLSuperClustersTableProducer",
    skipNonExistingSrc=cms.bool(True),
    src=cms.InputTag("particleFlowSuperClusterHGCal"),  # Offline collection (default)
    cut=cms.string(""),
    name=cms.string("TICLSuperClusters"),
    doc=cms.string("Offline TICL SuperClusters"),
    singleton=cms.bool(False),  # the number of entries is variable
    variables=cms.PSet(
        raw_energy=Var("rawEnergy", "float",
                       doc="Raw Energy of the SuperCluster [GeV]"),
        energy=Var("energy", "float",
                   doc="Regressed Energy of SuperCluster [GeV]"),
        corrected_energy=Var(
            "correctedEnergy", "float", doc="Corrected energy of the SuperCluster [GeV]"),
        position_x=Var("position.x", "float",
                       doc="SuperCluster position x [cm]"),
        position_y=Var("position.y", "float",
                       doc="SuperCluster position y [cm]"),
        position_z=Var("position.z", "float",
                       doc="SuperCluster barycenter z [cm]"),
        position_eta=Var("position.eta", "float",
                         doc="SuperCluster position pseudorapidity"),
        position_phi=Var("position.phi", "float",
                         doc="SuperCluster position phi"),
    ),
)
