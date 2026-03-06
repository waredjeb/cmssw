import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import *
from RecoHGCal.Configuration.hgcalTracksters_cfi import *
from RecoHGCal.Configuration.hgcalTICLCandidates_cfi import *
from RecoHGCal.Configuration.hgcalTICLSuperClusters_cfi import *
from RecoHGCal.Configuration.hgcalLayerClusters_cfi import *
from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5

######################################
# Offline HGCAL NanoAOD Tables
######################################

# Store HGCAL reconstruction objects
OfflineHGCalTables = cms.Sequence(
    hgcalTrackstersTableSequence
    + ticlCandidateTable
    + ticlCandidateExtraTable
)

# Add ticlSuperClustersTable only with ticl_v5 modifier
ticl_v5.toReplaceWith(
    OfflineHGCalTables,
    OfflineHGCalTables.copy() + ticlSuperClustersTable
)

# Store additional validation objects (SimTracksters, LayerClusters, associations)
OfflineHGCalValidationTables = cms.Sequence(
    hgcalTiclAssociationsTableSequence
    + hgcalSimTracksterSequence
    + ticlSimCandidateTable
    + ticlSimCandidateExtraTable
    + hgcalLayerClustersTableSequence
)

######################################
# Sequences for different NanoAOD flavours
######################################

# Offline HGCAL NanoAOD (NANO:@HGCAL) - reconstruction objects only
hgcalNanoSequence = cms.Sequence(
    OfflineHGCalTables
)

# Offline HGCAL NanoAOD with validation info (NANO:@HGCALVal) - includes MC/validation objects
hgcalNanoValidationSequence = cms.Sequence(
    OfflineHGCalTables
    + OfflineHGCalValidationTables
)

######################################
# Customization
######################################

def hgcalNanoCustomize(process):
    """
    Customization function for offline HGCAL NanoAOD.
    This function is called when producing NanoAOD with HGCAL content.
    """
    if hasattr(process, "NANOAODSIMoutput"):
        # Keep all HGCAL NanoAOD flat tables
        process.NANOAODSIMoutput.outputCommands.append(
            "keep nanoaodFlatTable_*Table*_*_*"
        )

    if hasattr(process, "NANOAODoutput"):
        # Keep all HGCAL NanoAOD flat tables
        process.NANOAODoutput.outputCommands.append(
            "keep nanoaodFlatTable_*Table*_*_*"
        )

    return process
