import FWCore.ParameterSet.Config as cms
from RecoHGCal.Configuration.hgcalTracksters_cfi import createTracksterTables
from Validation.HGCalValidation.HLT_TICLIterLabels_cff import hltTiclIterLabels

hltSimTrackstersLabels = [
    'hltTiclSimTracksters', 'hltTiclSimTrackstersfromCPs']

# Create HLT trackster tables 
_hltProducers = createTracksterTables(hltTiclIterLabels, hltSimTrackstersLabels, collectionPrefix="hlt")

for name, producer in _hltProducers.items():
    globals()[name] = producer.clone()

from RecoHGCal.Configuration.hgcalTracksters_cfi import SimCl2CPOneToOneFlatTable
hltSimCl2CPOneToOneFlatTable = SimCl2CPOneToOneFlatTable.clone()

tracksterTableProducers = []
hltTrackstersAssociationOneToManyTableProducers = []
simTracksterTableProducers = []

for name, producer in _hltProducers.items():
    if "AssociationTableProducer" in name and "SimCl2CP" not in name:
        hltTrackstersAssociationOneToManyTableProducers.append(globals()[name])
    elif "SimTrackster" in name or "fromCPs" in name:
        simTracksterTableProducers.append(globals()[name])
    elif "TableProducer" in name and "Association" not in name and "SimCl2CP" not in name:
        tracksterTableProducers.append(globals()[name])

hltTrackstersTableSequence = cms.Sequence(sum(tracksterTableProducers, cms.Sequence()))
hltTiclAssociationsTableSequence = cms.Sequence(sum(hltTrackstersAssociationOneToManyTableProducers, cms.Sequence()))
hltSimTracksterSequence = cms.Sequence(sum(simTracksterTableProducers, cms.Sequence()))

hltTiclAssociationsTableSequence += hltSimCl2CPOneToOneFlatTable
