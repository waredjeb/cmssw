import FWCore.ParameterSet.Config as cms
from RecoHGCal.Configuration.hgcalTracksters_cfi import createTracksterTables
from Validation.HGCalValidation.HLT_TICLIterLabels_cff import hltTiclIterLabels

hltSimTrackstersLabels = [
    'hltTiclSimTracksters', 'hltTiclSimTrackstersfromCPs']

# Create HLT trackster tables using factory function from RecoHGCal
_hltProducers = createTracksterTables(hltTiclIterLabels, hltSimTrackstersLabels, collectionPrefix="hlt")

# Assign all producers to module globals
for name, producer in _hltProducers.items():
    globals()[name] = producer

# Build sequences for organizing producers
tracksterTableProducers = []
hltTrackstersAssociationOneToManyTableProducers = []
simTracksterTableProducers = []

for name, producer in _hltProducers.items():
    if "AssociationTableProducer" in name and "SimCl2CP" not in name:
        hltTrackstersAssociationOneToManyTableProducers.append(producer)
    elif "SimTrackster" in name or "fromCPs" in name:
        simTracksterTableProducers.append(producer)
    elif "TableProducer" in name and "Association" not in name and "SimCl2CP" not in name:
        tracksterTableProducers.append(producer)

# Create sequences
hltTrackstersTableSequence = cms.Sequence(sum(tracksterTableProducers, cms.Sequence()))
hltTiclAssociationsTableSequence = cms.Sequence(sum(hltTrackstersAssociationOneToManyTableProducers, cms.Sequence()))
hltSimTracksterSequence = cms.Sequence(sum(simTracksterTableProducers, cms.Sequence()))

# Add SimCl2CP producer
hltSimCl2CPOneToOneFlatTable = _hltProducers['hltSimCl2CPOneToOneFlatTable']
hltTiclAssociationsTableSequence += hltSimCl2CPOneToOneFlatTable
