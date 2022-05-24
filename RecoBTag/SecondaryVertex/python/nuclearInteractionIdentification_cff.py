import FWCore.ParameterSet.Config as cms

from RecoBTag.SecondaryVertex.nuclearInteractionIdentifier_cfi import *
from RecoBTag.SecondaryVertex.vertexAndTracksCleaner_cfi import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *


#===========================================
# NI rejection and tracks cleaning procedure
#
# short descriptions of the different versions (0 to 7)
# version 0: identify solely based on position
# version 1: identify based on position broader position then geometry
# version 2-3: identify based on position & Nctau (=flightDistance2D/(gamma*Bctau) with gamma = pt/mass and Bctau = 0.05 cm
# version 4-5: identify based on position & ntracks 
# version 4-7: same as versions 0-3, but first run IVF with relaxed cuts (see below)
#
#===========================================




# IVF run with relaxed cuts to identify more NIs

inclusiveVertexFinderRelaxed = inclusiveCandidateVertexFinder.clone(
    minHits = (0),
    maximumLongitudinalImpactParameter = -9999.9, # 3.0 by Dominik
    vertexMinAngleCosine = -9999.9, # 0.7 by Dominik
    vertexMinDLen2DSig = -9999.9,
    vertexMinDLenSig = -9999.9
)

vertexMergerRelaxed = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderRelaxed") )

trackVertexArbitratorRelaxed = candidateVertexArbitrator.clone( secondaryVertices = cms.InputTag("vertexMergerRelaxed"))

inclusiveSecondaryVerticesRelaxed = candidateVertexMerger.clone(
    secondaryVertices = cms.InputTag("trackVertexArbitratorRelaxed")
,    maxFraction = cms.double(9999.9), # 0.2 by Dominik
    minSignificance = cms.double(-9999.9) # 10. by Dominik
)

inclusiveCandidateVertexingRelaxed = cms.Sequence(inclusiveVertexFinderRelaxed*vertexMergerRelaxed*trackVertexArbitratorRelaxed*inclusiveSecondaryVerticesRelaxed)

# NI identifiers

nuclearInteractionIdentifier0 = nuclearInteractionCandIdentifier.clone(
    selection = cms.PSet(
        nuclearInteractionCandIdentifier.selection,
#        position = cms.vdouble(2.65, 3.22, 3.52, 5.11, 6.64, 8.01, 9.53, 10.64) # old geometry
        position = cms.vdouble(2.17, 2.25, 2.694, 3.481, 6.597, 7.278, 10.700, 11.344, 15.797, 16.430)
#        position = cms.vdouble(1.9, 2.5, 2.694, 3.481, 6.597, 7.278, 10.700, 11.344, 15.797, 16.430) # broader beam pipe for checks with just removing beam pipe
#        position = cms.vdouble(5.0, 16.0) # forward region
    )
)

nuclearInteractionIdentifier1 = nuclearInteractionCandIdentifier.clone(
    selection = cms.PSet(
# best cut for deepCSV
        nuclearInteractionIdentifier0.selection,
        minNtracks = cms.int32(0),
        maxNtracks = cms.int32(2),
        # minNctau = cms.double(3.5),
        # minMass = cms.double(0.),
        # maxMass = cms.double(1.5),
        # distToNI = cms.double(0.1)
    )
)

nuclearInteractionIdentifier2 = nuclearInteractionCandIdentifier.clone(
    selection = cms.PSet(
        nuclearInteractionIdentifier0.selection,
#        minNtracks = cms.int32(0),
#        maxNtracks = cms.int32(2),
        minNctau = cms.double(3.5)
    )
)


nuclearInteractionIdentifier3 = nuclearInteractionCandIdentifier.clone(
    selection = cms.PSet(
        nuclearInteractionIdentifier0.selection,
        # minNtracks = cms.int32(0),
        # maxNtracks = cms.int32(2),
        minMass = cms.double(0.),
        maxMass = cms.double(1.5)
        )
)

nuclearInteractionIdentifier4 = nuclearInteractionCandIdentifier.clone(
    selection = cms.PSet(
        nuclearInteractionIdentifier0.selection,
        # minNtracks = cms.int32(0),
        # maxNtracks = cms.int32(2),
        distToNI = cms.double(0.1)
    )
)

nuclearInteractionIdentifier5 = nuclearInteractionCandIdentifier.clone(
    selection = cms.PSet(
        nuclearInteractionIdentifier0.selection,
        minNtracks = cms.int32(0),
        maxNtracks = cms.int32(3),
        # minNctau = cms.double(3.5),
        # minMass = cms.double(0.),
        # maxMass = cms.double(1.5)
    )
)

nuclearInteractionIdentifier6 = nuclearInteractionCandIdentifier.clone(
    selection = cms.PSet(
        nuclearInteractionIdentifier0.selection,
        minNtracks = cms.int32(0),
        maxNtracks = cms.int32(4),
        # minNctau = cms.double(3.5),
        # distToNI = cms.double(0.1)
    )
)

nuclearInteractionIdentifier7 = nuclearInteractionCandIdentifier.clone(
    selection = cms.PSet(
        nuclearInteractionIdentifier0.selection,
        # minNtracks = cms.int32(0),
        # maxNtracks = cms.int32(2),
        distToNI = cms.double(0.1),
        minMass = cms.double(0.),
        maxMass = cms.double(3.5)
        )
)

# vertex and pfcandidates cleaning steps

vertexAndTracksCleaned0 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier0")
vertexAndTracksCleaned1 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier1")
vertexAndTracksCleaned2 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier2")
vertexAndTracksCleaned3 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier3")
vertexAndTracksCleaned4 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier4")
vertexAndTracksCleaned5 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier5")
vertexAndTracksCleaned6 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier6")
vertexAndTracksCleaned7 = vertexAndTracksCandCleaned.clone(veto = "nuclearInteractionIdentifier7")

# re-run IVF

inclusiveVertexFinderCleaned0 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned0"))
inclusiveVertexFinderCleaned1 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned1"))
inclusiveVertexFinderCleaned2 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned2"))
inclusiveVertexFinderCleaned3 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned3"))
inclusiveVertexFinderCleaned4 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned4"))
inclusiveVertexFinderCleaned5 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned5"))
inclusiveVertexFinderCleaned6 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned6"))
inclusiveVertexFinderCleaned7 = inclusiveCandidateVertexFinder.clone(tracks = cms.InputTag("vertexAndTracksCleaned7"))

vertexMergerCleaned0 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned0"))
vertexMergerCleaned1 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned1"))
vertexMergerCleaned2 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned2"))
vertexMergerCleaned3 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned3"))
vertexMergerCleaned4 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned4"))
vertexMergerCleaned5 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned5"))
vertexMergerCleaned6 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned6"))
vertexMergerCleaned7 = candidateVertexMerger.clone( secondaryVertices = cms.InputTag("inclusiveVertexFinderCleaned7"))

trackVertexArbitratorCleaned0 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned0"), secondaryVertices = cms.InputTag("vertexMergerCleaned0"))
trackVertexArbitratorCleaned1 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned1"), secondaryVertices = cms.InputTag("vertexMergerCleaned1"))
trackVertexArbitratorCleaned2 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned2"), secondaryVertices = cms.InputTag("vertexMergerCleaned2"))
trackVertexArbitratorCleaned3 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned3"), secondaryVertices = cms.InputTag("vertexMergerCleaned3"))
trackVertexArbitratorCleaned4 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned4"), secondaryVertices = cms.InputTag("vertexMergerCleaned4"))
trackVertexArbitratorCleaned5 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned5"), secondaryVertices = cms.InputTag("vertexMergerCleaned5"))
trackVertexArbitratorCleaned6 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned6"), secondaryVertices = cms.InputTag("vertexMergerCleaned6"))
trackVertexArbitratorCleaned7 = candidateVertexArbitrator.clone(tracks = cms.InputTag("vertexAndTracksCleaned7"), secondaryVertices = cms.InputTag("vertexMergerCleaned7"))

inclusiveSecondaryVerticesCleaned0 = candidateVertexMerger.clone(
    secondaryVertices = cms.InputTag("trackVertexArbitratorCleaned0"),
    maxFraction = cms.double(0.2),
    minSignificance = cms.double(10.)
)

inclusiveSecondaryVerticesCleaned1 = inclusiveSecondaryVerticesCleaned0.clone(secondaryVertices = cms.InputTag("trackVertexArbitratorCleaned1"))
inclusiveSecondaryVerticesCleaned2 = inclusiveSecondaryVerticesCleaned0.clone(secondaryVertices = cms.InputTag("trackVertexArbitratorCleaned2"))
inclusiveSecondaryVerticesCleaned3 = inclusiveSecondaryVerticesCleaned0.clone(secondaryVertices = cms.InputTag("trackVertexArbitratorCleaned3"))
inclusiveSecondaryVerticesCleaned4 = inclusiveSecondaryVerticesCleaned0.clone(secondaryVertices = cms.InputTag("trackVertexArbitratorCleaned4"))
inclusiveSecondaryVerticesCleaned5 = inclusiveSecondaryVerticesCleaned0.clone(secondaryVertices = cms.InputTag("trackVertexArbitratorCleaned5"))
inclusiveSecondaryVerticesCleaned6 = inclusiveSecondaryVerticesCleaned0.clone(secondaryVertices = cms.InputTag("trackVertexArbitratorCleaned6"))
inclusiveSecondaryVerticesCleaned7 = inclusiveSecondaryVerticesCleaned0.clone(secondaryVertices = cms.InputTag("trackVertexArbitratorCleaned7"))

# all NI rejection sequences

nuclearInteractionsRemoved0 = cms.Sequence(
    inclusiveCandidateVertexing *
    nuclearInteractionIdentifier0 *
    vertexAndTracksCleaned0 *
    inclusiveVertexFinderCleaned0 *
    vertexMergerCleaned0 *
    trackVertexArbitratorCleaned0 *
    inclusiveSecondaryVerticesCleaned0
)


nuclearInteractionsRemoved1 = cms.Sequence(
    inclusiveCandidateVertexing *
    nuclearInteractionIdentifier1 *
    vertexAndTracksCleaned1 *
    inclusiveVertexFinderCleaned1 *
    vertexMergerCleaned1 *
    trackVertexArbitratorCleaned1 *
    inclusiveSecondaryVerticesCleaned1
)

nuclearInteractionsRemoved2 = cms.Sequence(
    inclusiveCandidateVertexing *
    nuclearInteractionIdentifier2 *
    vertexAndTracksCleaned2 *
    inclusiveVertexFinderCleaned2 *
    vertexMergerCleaned2 *
    trackVertexArbitratorCleaned2 *
    inclusiveSecondaryVerticesCleaned2
)

nuclearInteractionsRemoved3 = cms.Sequence(
    inclusiveCandidateVertexing *
    nuclearInteractionIdentifier3 *
    vertexAndTracksCleaned3 *
    inclusiveVertexFinderCleaned3 *
    vertexMergerCleaned3 *
    trackVertexArbitratorCleaned3 *
    inclusiveSecondaryVerticesCleaned3
)

nuclearInteractionsRemoved4 = cms.Sequence(
    inclusiveCandidateVertexingRelaxed *
    nuclearInteractionIdentifier4 *
    vertexAndTracksCleaned4 *
    inclusiveVertexFinderCleaned4 *
    vertexMergerCleaned4 *
    trackVertexArbitratorCleaned4 *
    inclusiveSecondaryVerticesCleaned4
)

nuclearInteractionsRemoved5 = cms.Sequence(
    inclusiveCandidateVertexingRelaxed *
    nuclearInteractionIdentifier5 *
    vertexAndTracksCleaned5 *
    inclusiveVertexFinderCleaned5 *
    vertexMergerCleaned5 *
    trackVertexArbitratorCleaned5 *
    inclusiveSecondaryVerticesCleaned5
)

nuclearInteractionsRemoved6 = cms.Sequence(
    inclusiveCandidateVertexingRelaxed *
    nuclearInteractionIdentifier6 *
    vertexAndTracksCleaned6 *
    inclusiveVertexFinderCleaned6 *
    vertexMergerCleaned6 *
    trackVertexArbitratorCleaned6 *
    inclusiveSecondaryVerticesCleaned6
)

nuclearInteractionsRemoved7 = cms.Sequence(
    inclusiveCandidateVertexingRelaxed *
    nuclearInteractionIdentifier7 *
    vertexAndTracksCleaned7 *
    inclusiveVertexFinderCleaned7 *
    vertexMergerCleaned7 *
    trackVertexArbitratorCleaned7 *
    inclusiveSecondaryVerticesCleaned7
)
