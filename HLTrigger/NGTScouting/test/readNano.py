import uproot
import awkward as ak
import numpy as np
import argparse
import sys

def main():
    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("--validation", action="store_true", help="Check Sim branches")
    parser.add_argument("--ticlv", choices=["v4", "v5"], default="v5", help="TICL version")
    args = parser.parse_args()

    file_path = args.filename
    ticlVersion = args.ticlv
    validation = args.validation
    tree_name = "Events"

    # Open the ROOT file and load the TTree
    try:
        with uproot.open(file_path) as file:
            tree = file[tree_name]
            print(tree.keys())

            events = tree.arrays(library="ak")
            print("Fields:", "\n\t".join(events[0].fields))
    except Exception as e:
        print(f"Error opening file or reading tree: {e}")
        sys.exit(1)

    # Event loop
    for i, event in enumerate(events):
        print(f"Processing event {i}")
        ## CLUE3D Tracksters ##
        print("Found {} CLUE3D Tracksters".format(event.nhltTiclTrackstersCLUE3DHigh))
        for t_idx in range(event.nhltTiclTrackstersCLUE3DHigh):
            offset = event.hltTiclTrackstersCLUE3DHigh_ohltTiclTrackstersCLUE3DHighvertices[t_idx]
            count = event.hltTiclTrackstersCLUE3DHigh_nhltTiclTrackstersCLUE3DHighvertices[t_idx]
            vertices = event.hltTiclTrackstersCLUE3DHighvertices_vertices[offset : offset + count]
            vertex_multiplicity = event.hltTiclTrackstersCLUE3DHighvertices_vertex_mult[offset : offset + count]
            print(
                t_idx,
                list(zip(vertices, vertex_multiplicity)),
                event.hltTiclTrackstersCLUE3DHigh_raw_energy[t_idx],
            )

        print("Exploring connections, scores, and sharedEnergy")
        print("Connections for {} objects".format(event.nSimTSCP2hltTiclTrackstersCLUE3DHighByHitsLinks))
        try:  # Offset pattern
            offset = 0
            for obj_idx in range(event.nSimTSCP2hltTiclTrackstersCLUE3DHighMergeByHits - 1):
                next_offset = event.SimTSCP2hltTiclTrackstersCLUE3DHighMergeByHits_oSimTSCP2hltTiclTrackstersCLUE3DHighByHitsLinks[
                    obj_idx + 1
                ]
                elements = event.SimTSCP2hltTiclTrackstersCLUE3DHighHitsLinks_index[offset:next_offset]
                scores = event.SimTSCP2hltTiclTrackstersCLUE3DHighHitsLinks_score[offset:next_offset]
                sharedEnergy = event.SimTSCP2hltTiclTrackstersCLUE3DHighHitsLinks_sharedEnergy[
                    offset:next_offset
                ]
                if len(elements) > 0:
                    print("Offset ", obj_idx, elements, scores, sharedEnergy)
                offset = next_offset
        except AttributeError as e:
            print(f"An AttributeError occurred (Offset): {e}")

        try:  # Count pattern
            offset = 0
            for obj_idx in range(event.nSimTSCP2hltTiclTrackstersCLUE3DHighMergeByHits):
                count = event.SimTSCP2hltTiclTrackstersCLUE3DHighMergeByHits_nSimTSCP2hltTiclTrackstersCLUE3DHighByHitsLinks[obj_idx]
                elements = event.SimTSCP2hltTiclTrackstersCLUE3DHighByHitsLinks_index[offset : offset + count]
                scores = event.SimTSCP2hltTiclTrackstersCLUE3DHighByHitsLinks_score[offset : offset + count]
                sharedEnergy = event.SimTSCP2hltTiclTrackstersCLUE3DHighByHitsLinks_sharedEnergy[
                    offset : offset + count
                ]
                if len(elements) > 0:
                    print("Count ", obj_idx, elements, scores, sharedEnergy)
                offset += count
        except AttributeError as e:
            print(f"An AttributeError occurred (Count): {e}")

        ## TICLv5 Collections ##
        if(ticlVersion == "v5"):
            print("Found {} TICLCandidate tracksters".format(event.nhltTiclCandidate))
            for t_idx in range(event.nhltTiclCandidate):
                offset = event.hltTiclCandidate_ohltTiclCandidatevertices[t_idx]
                count = event.hltTiclCandidate_nhltTiclCandidatevertices[t_idx]
                vertices = event.hltTiclCandidatevertices_vertices[offset : offset + count]
                vertex_multiplicity = event.hltTiclCandidatevertices_vertex_mult[offset : offset + count]
                print(
                    t_idx,
                    list(zip(vertices, vertex_multiplicity)),
                    event.hltTiclCandidate_raw_energy[t_idx],
                )
    
            print("Exploring connections, scores, and sharedEnergy")
            print("Connections for {} objects".format(event.nSimTSCP2hltTiclCandidateByHitsLinks))
            try:  # Offset pattern
                offset = 0
                for obj_idx in range(event.nSimTSCP2hltTiclCandidateMergeByHits - 1):
                    next_offset = event.SimTSCP2hltTiclCandidateMergeByHits_oSimTSCP2hltTiclCandidateHitsLinks[
                        obj_idx + 1
                    ]
                    elements = event.SimTSCP2hltTiclCandidateHitsLinks_index[offset:next_offset]
                    scores = event.SimTSCP2hltTiclCandidateHitsLinks_score[offset:next_offset]
                    sharedEnergy = event.SimTSCP2hltTiclCandidateHitsLinks_sharedEnergy[
                        offset:next_offset
                    ]
                    if len(elements) > 0:
                        print("Offset ", obj_idx, elements, scores, sharedEnergy)
                    offset = next_offset
            except AttributeError as e:
                print(f"An AttributeError occurred (Offset): {e}")
    
            try:  # Count pattern
                offset = 0
                for obj_idx in range(event.nSimTSCP2hltTiclCandidateMergeByHits):
                    count = event.SimTSCP2hltTiclTrackstersCLUE3DHighMergeByHits_nSimTSCP2hltTiclTrackstersCLUE3DHighByHitsLinks[obj_idx]
                    elements = event.SimTSSC2TShltTiclCandidateByHitsLinks_index[offset : offset + count]
                    scores = event.SimTSSC2TShltTiclCandidateByHitsLinks_score[offset : offset + count]
                    sharedEnergy = event.SimTSSC2TShltTiclCandidateByHitsLinks_sharedEnergy[
                        offset : offset + count
                    ]
                    if len(elements) > 0:
                        print("Count ", obj_idx, elements, scores, sharedEnergy)
                    offset += count
            except AttributeError as e:
                print(f"An AttributeError occurred (Count): {e}")

            if(validation):
                print("Found {} simTICLCandidates".format(event.nhltSimTICLCandidates))
                for sim_idx in range(event.nhltSimTICLCandidates):
                    trackIdx = event.hltSimTICLCandidates_trackIdx[sim_idx]
                    if(trackIdx >= 0):
                        track_pt = event.hltGeneralTrack_pt[trackIdx]
                    else:
                        track_pt = np.nan;
                    print(
                        sim_idx,
                        trackIdx,
                        track_pt,
                        event.hltSimTICLCandidates_raw_energy[sim_idx],
                    )

        try:
            for i in range(event.nSimCl2CPWithFraction):
                print(
                    "SimCl {} is linked to CP {} with fraction {}".format(
                        i,
                        event.SimCl2CPWithFraction_index[i],
                        event.SimCl2CPWithFraction_fraction[i],
                    )
                )
        except AttributeError as e:
            print(f"An AttributeError occurred (SimCl2CPWithFraction): {e}")

if __name__ == "__main__":
    main()
