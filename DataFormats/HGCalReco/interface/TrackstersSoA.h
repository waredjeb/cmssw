#ifndef DataFormats_HGCalReco_interface_TrackstersSoA_h
#define DataFormats_HGCalReco_interface_TrackstersSoA_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoABlocks.h"

GENERATE_SOA_LAYOUT(GNNNodeSoALayout,
                    // columns: one value per element
                    SOA_COLUMN(float, barycenter_x),
                    SOA_COLUMN(float, barycenter_y),
                    SOA_COLUMN(float, barycenter_z),
                    SOA_COLUMN(float, barycenter_eta),
                    SOA_COLUMN(float, barycenter_phi),
                    SOA_COLUMN(float, eigenvector0_x),
                    SOA_COLUMN(float, eigenvector0_y),
                    SOA_COLUMN(float, eigenvector0_z),
                    SOA_COLUMN(float, eigenvalue1),
                    SOA_COLUMN(float, eigenvalue2),
                    SOA_COLUMN(float, eigenvalue3),
                    SOA_COLUMN(float, sigmasPCA1),
                    SOA_COLUMN(float, sigmasPCA2),
                    SOA_COLUMN(float, sigmasPCA3),
                    SOA_COLUMN(float, num_LCs),
                    SOA_COLUMN(float, num_hits),
                    SOA_COLUMN(float, raw_energy),
                    SOA_COLUMN(float, raw_em_energy),
                    SOA_COLUMN(float, photon_prob),
                    SOA_COLUMN(float, electron_prob),
                    SOA_COLUMN(float, muon_prob),
                    SOA_COLUMN(float, neutral_pion_prob),
                    SOA_COLUMN(float, charged_hadron_prob),
                    SOA_COLUMN(float, neutral_hadron_prob),
                    SOA_COLUMN(float, z_min),
                    SOA_COLUMN(float, z_max),
                    SOA_COLUMN(float, LC_density),
                    SOA_COLUMN(float, trackster_density),
                    SOA_COLUMN(float, time));

GENERATE_SOA_LAYOUT(GNNEdgeSoALayout,
                    SOA_COLUMN(float, max_raw_energy),
                    SOA_COLUMN(float, raw_energy),
                    SOA_COLUMN(float, barycenter_z),
                    SOA_COLUMN(float, barycenter_xy),
                    SOA_COLUMN(float, eigenvector0),
                    SOA_COLUMN(float, time));

GENERATE_SOA_LAYOUT(GNNEdgeIndexSoALayout, SOA_COLUMN(long, in), SOA_COLUMN(long, out));

// SoA-by-blocks layout combining the node, edge and edge-index blocks into a single
// portable product (replaces the removed PortableMultiCollection<node, edge, edgeIndex>).
GENERATE_SOA_BLOCKS(GNNTrackstersBlocksLayout,
                    SOA_BLOCK(nodes, GNNNodeSoALayout),
                    SOA_BLOCK(edges, GNNEdgeSoALayout),
                    SOA_BLOCK(edgeIndex, GNNEdgeIndexSoALayout))

GENERATE_SOA_LAYOUT(GNNOutputSoALayout, SOA_COLUMN(float, score));

GENERATE_SOA_LAYOUT(GNNPostprocessingSoALayout,
                    SOA_COLUMN(float, score),
                    SOA_COLUMN(long, in),
                    SOA_COLUMN(long, out),
                    SOA_COLUMN(float, max_raw_energy));

using GNNNodeSoA = GNNNodeSoALayout<>;
using GNNEdgeSoA = GNNEdgeSoALayout<>;
using GNNEdgeIndexSoA = GNNEdgeIndexSoALayout<>;
using GNNTrackstersBlocks = GNNTrackstersBlocksLayout<>;
using GNNOutputSoA = GNNOutputSoALayout<>;
using GNNPostprocessingSoA = GNNPostprocessingSoALayout<>;

#endif
