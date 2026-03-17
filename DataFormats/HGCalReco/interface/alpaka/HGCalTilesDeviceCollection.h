#ifndef DataFormats_HGCalReco_interface_alpaka_HGCalTilesDeviceCollection_h
#define DataFormats_HGCalReco_interface_alpaka_HGCalTilesDeviceCollection_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

// Flattened tiles structure for device-friendly access
// Original CPU tiles: vector<vector<vector<int>>> (layer -> etaphi bin -> cluster indices)
// Flattened: offsets array + content array
GENERATE_SOA_LAYOUT(HGCalTilesSoALayout,
                    // Tile structure
                    SOA_COLUMN(int, tileOffsets),  // prefix sum, size = nLayers * nEtaBins * nPhiBins + 1
                    SOA_COLUMN(int, tileContent),  // flattened cluster indices

                    // Tile dimensions
                    SOA_SCALAR(int, nLayers),
                    SOA_SCALAR(int, nEtaBins),
                    SOA_SCALAR(int, nPhiBins),

                    // Binning parameters (for bin calculation on device)
                    SOA_SCALAR(float, etaMin),
                    SOA_SCALAR(float, etaMax),
                    SOA_SCALAR(float, phiMin),
                    SOA_SCALAR(float, phiMax)
)

using HGCalTilesSoA = HGCalTilesSoALayout<>;

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using HGCalTilesDeviceCollection = PortableCollection<HGCalTilesSoA>;
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // DataFormats_HGCalReco_interface_alpaka_HGCalTilesDeviceCollection_h
