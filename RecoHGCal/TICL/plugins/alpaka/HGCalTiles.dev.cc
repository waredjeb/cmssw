#include "HGCalTiles.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace ticl {

    // Kernel to count clusters per tile
    struct CountClustersPerTileKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(const TAcc& acc,
                                   CLUE3DStateSoA::ConstView clusters,
                                   HGCalTilesView tilesView,
                                   int nClusters) const {
        for (int i : cms::alpakatools::uniform_elements(acc, nClusters)) {
          float eta = clusters.eta(i);
          float phi = clusters.phi(i);
          int layer = clusters.layer(i);

          int globalBin = tilesView.getGlobalBin(layer, eta, phi);
          alpaka::atomicAdd(acc, &tilesView.offsets[globalBin], 1, alpaka::hierarchy::Blocks{});
        }
      }
    };

    // Kernel to perform exclusive scan (prefix sum)
    struct ExclusiveScanKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(const TAcc& acc, int* data, int size) const {
        // Simple serial scan for now (could be optimized with parallel scan)
        if (alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0] == 0) {
          int sum = 0;
          for (int i = 0; i < size; ++i) {
            int temp = data[i];
            data[i] = sum;
            sum += temp;
          }
          data[size] = sum;  // Store total
        }
      }
    };

    // Kernel to fill tile indexes
    struct FillTileIndexesKernel {
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(const TAcc& acc,
                                   CLUE3DStateSoA::ConstView clusters,
                                   HGCalTilesView tilesView,
                                   int* binCounters,  // temporary counter per bin
                                   int nClusters) const {
        for (int i : cms::alpakatools::uniform_elements(acc, nClusters)) {
          float eta = clusters.eta(i);
          float phi = clusters.phi(i);
          int layer = clusters.layer(i);

          int globalBin = tilesView.getGlobalBin(layer, eta, phi);

          // Atomically get position and increment
          int pos = alpaka::atomicAdd(acc, &binCounters[globalBin], 1, alpaka::hierarchy::Blocks{});
          int offset = tilesView.offsets[globalBin];

          // Write cluster index
          tilesView.indexes[offset + pos] = i;
        }
      }
    };

    template <typename TQueue>
    void HGCalTiles::fill(TQueue& queue,
                         const CLUE3DStateDeviceCollection& clusters,
                         int nClusters) {
      // Step 1: Zero out offsets
      alpaka::memset(queue, offsets_, 0);

      // Step 2: Count clusters per tile (results go into offsets temporarily)
      auto workDiv1 = cms::alpakatools::make_workdiv<Acc1D>(nClusters, 256);
      alpaka::exec<Acc1D>(queue, workDiv1, CountClustersPerTileKernel{}, clusters.const_view(), view_, nClusters);

      // Step 3: Exclusive scan to get offsets
      auto workDiv2 = cms::alpakatools::make_workdiv<Acc1D>(1, 1);
      alpaka::exec<Acc1D>(
          queue, workDiv2, ExclusiveScanKernel{}, alpaka::getPtrNative(offsets_), nTiles_ + 1);

      // Step 4: Allocate temporary counter buffer (initialized to 0)
      auto binCounters = cms::alpakatools::make_device_buffer<int[]>(queue, nTiles_);
      alpaka::memset(queue, binCounters, 0);

      // Step 5: Fill indexes array
      auto workDiv3 = cms::alpakatools::make_workdiv<Acc1D>(nClusters, 256);
      alpaka::exec<Acc1D>(
          queue, workDiv3, FillTileIndexesKernel{}, clusters.const_view(), view_, alpaka::getPtrNative(binCounters), nClusters);
    }

  }  // namespace ticl
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
