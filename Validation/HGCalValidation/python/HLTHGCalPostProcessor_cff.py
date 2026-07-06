import FWCore.ParameterSet.Config as cms

from Validation.HGCalValidation.PostProcessorHGCAL_cfi import postProcessorHGCALlayerclusters as _postProcessorHGCALlayerclusters
from Validation.HGCalValidation.PostProcessorHGCAL_cfi import postProcessorHGCALsimclusters as _postProcessorHGCALsimclusters
from Validation.HGCalValidation.PostProcessorHGCAL_cfi import postProcessorHGCALTracksters as _postProcessorHGCALTracksters
from Validation.HGCalValidation.PostProcessorHGCAL_cfi import postProcessorHGCALCandidates as _postProcessorHGCALCandidates 

from Validation.HGCalValidation.HLT_TICLIterLabels_cff import hltTiclIterLabelsPSet as _hltTiclIterLabelsPSet
from Validation.HGCalValidation.HLTHGCalValidator_cff import hltHgcalValidator as _hltHgcalValidator

hltPrefix = 'HLT/HGCAL/HGCalValidator/'
hltTracksterLabels = _hltTiclIterLabelsPSet.labels.copy()
hltTracksterLabels.extend(['hltTiclSimTracksters', 'hltTiclSimTracksters_fromCPs'])
# L1-seeded SimTracksters validated (byHits only) vs the unseeded SimTracksters:
# their efficiency/response plots live under the unseeded validator's dir too.
hltTracksterLabels.extend(['hltTiclSimTrackstersL1Seeded', 'hltTiclSimTrackstersL1Seeded_fromCPs'])

hltLcToCP_linking = _hltHgcalValidator.label_LCToCPLinking.value()
hltPostProcessorHGCALlayerclusters = _postProcessorHGCALlayerclusters.clone(
    subDirs = cms.untracked.vstring(hltPrefix + _hltHgcalValidator.label_layerClustersPlots.value() + '/' + hltLcToCP_linking),
)

hltSubdirsSim = [hltPrefix + _hltHgcalValidator.label_SimClusters.value() + '/'+iteration+'/' for iteration in hltTracksterLabels]
hltPostProcessorHGCALsimclusters = _postProcessorHGCALsimclusters.clone(
    subDirs = cms.untracked.vstring(hltSubdirsSim)
)

hltTSbyHits_CP = _hltHgcalValidator.label_TSbyHitsCP.value()
hltSubdirsTracksters = [hltPrefix+iteration+'/'+hltTSbyHits_CP for iteration in hltTracksterLabels]

hltTSbyLCs = _hltHgcalValidator.label_TSbyLCs.value()
hltSubdirsTracksters.extend(hltPrefix+iteration+'/'+hltTSbyLCs for iteration in hltTracksterLabels)

hltTSbyLCs_CP = _hltHgcalValidator.label_TSbyLCsCP.value()
hltSubdirsTracksters.extend(hltPrefix+iteration+'/'+hltTSbyLCs_CP for iteration in hltTracksterLabels)

hltTSbyHits = _hltHgcalValidator.label_TSbyHits.value()
hltSubdirsTracksters.extend(hltPrefix+iteration+'/'+hltTSbyHits for iteration in hltTracksterLabels)

hltPostProcessorHGCALTracksters = _postProcessorHGCALTracksters.clone(
    subDirs = cms.untracked.vstring(hltSubdirsTracksters)
)

hltNeutrals = ["photons", "neutral_pions", "neutral_hadrons"]
hltCharged = ["electrons", "muons", "charged_hadrons"]
hltSubDirsCandidates = [hltPrefix + _hltHgcalValidator.ticlCandidates.value() + "/" + c for cands in (hltNeutrals, hltCharged) for c in cands]

hltPostProcessorHGCALCandidates = _postProcessorHGCALCandidates.clone(
    subDirs = cms.untracked.vstring(hltSubDirsCandidates)
)

hltHcalValidatorPostProcessor = cms.Sequence(
    hltPostProcessorHGCALlayerclusters+
    hltPostProcessorHGCALsimclusters+
    hltPostProcessorHGCALTracksters+
    hltPostProcessorHGCALCandidates
)

# --- L1-seeded validator harvesting ---------------------------------------
# Same clients, pointed at the L1-seeded validator's dirName and trackster set
# (candidate plots are off in the L1-seeded validator, so no candidates client).
from Validation.HGCalValidation.HLTHGCalValidator_cff import hltHgcalValidatorL1Seeded as _hltHgcalValidatorL1Seeded
hltPrefixL1S = 'HLT/HGCAL/HGCalValidatorL1Seeded/'
hltTracksterLabelsL1S = [l for l in _hltTiclIterLabelsPSet.labels if l.endswith("L1Seeded")]
hltTracksterLabelsL1S.extend(['hltTiclSimTrackstersL1Seeded', 'hltTiclSimTrackstersL1Seeded_fromCPs'])

hltPostProcessorHGCALlayerclustersL1Seeded = _postProcessorHGCALlayerclusters.clone(
    subDirs = cms.untracked.vstring(hltPrefixL1S + _hltHgcalValidatorL1Seeded.label_layerClustersPlots.value() + '/' + hltLcToCP_linking),
)

hltSubdirsSimL1S = [hltPrefixL1S + _hltHgcalValidatorL1Seeded.label_SimClusters.value() + '/'+iteration+'/' for iteration in hltTracksterLabelsL1S]
hltPostProcessorHGCALsimclustersL1Seeded = _postProcessorHGCALsimclusters.clone(
    subDirs = cms.untracked.vstring(hltSubdirsSimL1S)
)

hltSubdirsTrackstersL1S = [hltPrefixL1S+iteration+'/'+hltTSbyHits_CP for iteration in hltTracksterLabelsL1S]
hltSubdirsTrackstersL1S.extend(hltPrefixL1S+iteration+'/'+hltTSbyLCs for iteration in hltTracksterLabelsL1S)
hltSubdirsTrackstersL1S.extend(hltPrefixL1S+iteration+'/'+hltTSbyLCs_CP for iteration in hltTracksterLabelsL1S)
hltSubdirsTrackstersL1S.extend(hltPrefixL1S+iteration+'/'+hltTSbyHits for iteration in hltTracksterLabelsL1S)

hltPostProcessorHGCALTrackstersL1Seeded = _postProcessorHGCALTracksters.clone(
    subDirs = cms.untracked.vstring(hltSubdirsTrackstersL1S)
)

hltHcalValidatorL1SeededPostProcessor = cms.Sequence(
    hltPostProcessorHGCALlayerclustersL1Seeded+
    hltPostProcessorHGCALsimclustersL1Seeded+
    hltPostProcessorHGCALTrackstersL1Seeded
)
