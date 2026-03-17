#ifndef DataFormats_HGCalReco_interface_CLUE3DStateSoA_h
#define DataFormats_HGCalReco_interface_CLUE3DStateSoA_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

// SoA layout for CLUE3D clustering algorithm state
// This stores per-cluster temporary state needed during pattern recognition
GENERATE_SOA_LAYOUT(CLUE3DStateSoALayout,
                    // Input cluster properties (copied/computed from layer clusters)
                    SOA_COLUMN(float, x),
                    SOA_COLUMN(float, y),
                    SOA_COLUMN(float, z),
                    SOA_COLUMN(float, eta),
                    SOA_COLUMN(float, phi),
                    SOA_COLUMN(float, r_over_absz),
                    SOA_COLUMN(float, radius),  // estimated cluster radius
                    SOA_COLUMN(float, energy),
                    SOA_COLUMN(int, cells),
                    SOA_COLUMN(int, layer),         // layer index
                    SOA_COLUMN(int, algoId),        // cluster algo ID
                    SOA_COLUMN(uint8_t, isSilicon),
                    SOA_COLUMN(int, layerClusterOriginalIdx),  // original index in input collection

                    // Density calculation results
                    SOA_COLUMN(float, rho),  // local density
                    SOA_COLUMN(float, z_extension),  // Z-span of density search window

                    // Distance to higher density results
                    SOA_COLUMN(float, delta_dist),   // distance to nearest higher (transverse)
                    SOA_COLUMN(int, delta_layer),    // distance in layers to nearest higher
                    SOA_COLUMN(int, nearestHigher_layer),  // layer of nearest higher density cluster
                    SOA_COLUMN(int, nearestHigher_idx),    // index in that layer of nearest higher

                    // Seed finding and assignment results
                    SOA_COLUMN(int, clusterIndex),  // assigned trackster index (-1 if not assigned)
                    SOA_COLUMN(uint8_t, isSeed),    // is this cluster a seed?
                    SOA_COLUMN(uint8_t, isOutlier)  // is this cluster an outlier?
)

using CLUE3DStateSoA = CLUE3DStateSoALayout<>;

#endif  // DataFormats_HGCalReco_interface_CLUE3DStateSoA_h
