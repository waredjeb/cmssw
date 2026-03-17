#ifndef RecoHGCal_TICL_plugins_alpaka_HGCalTiles_h
#define RecoHGCal_TICL_plugins_alpaka_HGCalTiles_h

#include "DataFormats/HGCalReco/interface/CLUE3DStateSoA.h"
#include "DataFormats/HGCalReco/interface/alpaka/CLUE3DStateDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace ticl {

    // Lightweight view for device-side access
    struct HGCalTilesView {
      int* indexes;    // flattened cluster indices
      int* offsets;    // prefix sum for bin boundaries
      int nLayers;
      int nEtaBins;
      int nPhiBins;
      float etaMin;
      float etaMax;
      float phiMin;
      float phiMax;

      // Get bin index for given eta/phi/layer
      ALPAKA_FN_ACC inline int getGlobalBin(int layer, float eta, float phi) const {
        int etaBin = static_cast<int>((eta - etaMin) / (etaMax - etaMin) * nEtaBins);
        int phiBin = static_cast<int>((phi - phiMin) / (phiMax - phiMin) * nPhiBins);

        // Clamp to valid range
        if (etaBin < 0) etaBin = 0;
        if (etaBin >= nEtaBins) etaBin = nEtaBins - 1;
        if (phiBin < 0) phiBin = 0;
        if (phiBin >= nPhiBins) phiBin = nPhiBins - 1;

        return layer * nEtaBins * nPhiBins + etaBin * nPhiBins + phiBin;
      }

      // Get content of a tile
      ALPAKA_FN_ACC inline int* operator[](int globalBin) const {
        return indexes + offsets[globalBin];
      }

      // Get size of a tile
      ALPAKA_FN_ACC inline int size(int globalBin) const {
        return offsets[globalBin + 1] - offsets[globalBin];
      }
    };

    // High-level tiles manager
    class HGCalTiles {
    public:
      template <typename TQueue>
      HGCalTiles(TQueue const& queue, int nClusters, int nLayers, int nEtaBins, int nPhiBins)
          : indexes_(cms::alpakatools::make_device_buffer<int[]>(queue, nClusters)),
            offsets_(cms::alpakatools::make_device_buffer<int[]>(queue, nLayers * nEtaBins * nPhiBins + 1)),
            nLayers_(nLayers),
            nEtaBins_(nEtaBins),
            nPhiBins_(nPhiBins),
            nTiles_(nLayers * nEtaBins * nPhiBins) {

        // Initialize view
        view_.indexes = alpaka::getPtrNative(indexes_);
        view_.offsets = alpaka::getPtrNative(offsets_);
        view_.nLayers = nLayers;
        view_.nEtaBins = nEtaBins;
        view_.nPhiBins = nPhiBins;
        view_.etaMin = -3.5f;
        view_.etaMax = 3.5f;
        view_.phiMin = -M_PI;
        view_.phiMax = M_PI;
      }

      const HGCalTilesView& view() const { return view_; }

      // Fill tiles from cluster data
      template <typename TQueue>
      void fill(TQueue& queue,
                const CLUE3DStateDeviceCollection& clusters,
                int nClusters);

    private:
      cms::alpakatools::device_buffer<Device, int[]> indexes_;
      cms::alpakatools::device_buffer<Device, int[]> offsets_;
      int nLayers_;
      int nEtaBins_;
      int nPhiBins_;
      int nTiles_;
      HGCalTilesView view_;
    };

  }  // namespace ticl
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoHGCal_TICL_plugins_alpaka_HGCalTiles_h
