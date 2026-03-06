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

OfflineHGCalTables = cms.Sequence(
    hgcalTrackstersTableSequence
    + ticlCandidateTable
    + ticlCandidateExtraTable
)

# Store additional validation objects 
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

# Offline HGCAL NanoAOD with validation info (NANO:@HGCALVal) - includes sim objects and scores
hgcalNanoValidationSequence = cms.Sequence(
    OfflineHGCalTables
    + OfflineHGCalValidationTables
)

def hgcalNanoCustomize(process):
    """
    Customization function for offline HGCAL NanoAOD.
    This function is called when producing NanoAOD with HGCAL content.
    """
    if hasattr(process, "NANOAODSIMoutput"):
        process.NANOAODSIMoutput.outputCommands.append(
            "keep nanoaodFlatTable_*Table*_*_*"
        )

    if hasattr(process, "NANOAODoutput"):
        process.NANOAODoutput.outputCommands.append(
            "keep nanoaodFlatTable_*Table*_*_*"
        )

    return process
