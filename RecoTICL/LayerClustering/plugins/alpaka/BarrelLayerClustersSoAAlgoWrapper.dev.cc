#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "BarrelLayerClustersSoAAlgoWrapper.h"
#include "ConstantsForClusters.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;
  using namespace hgcal::constants;

  namespace {

    // Energy, number of hits, seed and energy-weighted position of each cluster.
    // Without sharing every rechit belongs to at most one cluster, with fraction 1.
    template <typename CAL>
    struct KernelAccumulate {
      ALPAKA_FN_ACC void operator()(const Acc1D& acc,
                                    reco::PFRecHitDeviceCollection::ConstView inputRecHits,
                                    HGCalSoARecHitsExtraDeviceCollection::ConstView inputClusters,
                                    reco::CaloClusterDeviceCollection::View outputs) const {
        for (auto i : cms::alpakatools::uniform_elements(acc, inputClusters.metadata().size())) {
          const auto clIdx = inputClusters[i].clusterIndex();
          if (clIdx == kInvalidCluster)
            continue;
          const auto energy = inputRecHits[i].energy();
          alpaka::atomicAdd(acc, &outputs.energy()[clIdx].energy(), energy);
          alpaka::atomicAdd(acc, &outputs.position()[clIdx].cells(), 1);
          alpaka::atomicAdd(acc, &outputs.position()[clIdx].x(), inputRecHits[i].x() * energy);
          alpaka::atomicAdd(acc, &outputs.position()[clIdx].y(), inputRecHits[i].y() * energy);
          alpaka::atomicAdd(acc, &outputs.position()[clIdx].z(), inputRecHits[i].z() * energy);
          if (inputClusters[i].isSeed() == 1) {
            outputs.indexes()[clIdx].seedID() = inputRecHits[i].detId();
            outputs.position()[clIdx].layer() = CAL::layer(inputRecHits, i);
          }
        }
      }
    };

    struct KernelFinalise {
      ALPAKA_FN_ACC void operator()(const Acc1D& acc,
                                    reco::CaloClusterDeviceCollection::View outputs,
                                    uint32_t numberOfClusters) const {
        for (auto clIdx : cms::alpakatools::uniform_elements(acc, numberOfClusters)) {
          // As in the legacy BarrelCLUEAlgo, a cluster without energy is placed at the origin
          const auto energy = outputs.energy()[clIdx].energy();
          const auto invEnergy = (energy != 0.f) ? 1.f / energy : 0.f;
          outputs.position()[clIdx].x() *= invEnergy;
          outputs.position()[clIdx].y() *= invEnergy;
          outputs.position()[clIdx].z() *= invEnergy;
          // Same defaults as reco::CaloCluster; the time is computed when converting to the legacy format
          outputs.energy()[clIdx].correctedEnergy() = -1.f;
          outputs.energy()[clIdx].correctedEnergyUncertainty() = -1.f;
          outputs.timing()[clIdx].time() = -99.f;
          outputs.timing()[clIdx].timeError() = -1.f;
          outputs.indexes()[clIdx].flags() = 0;
        }
      }
    };

  }  // namespace

  template <typename CAL>
  void BarrelLayerClustersSoAAlgoWrapper<CAL>::run(Queue& queue,
                                                   const unsigned int numberOfClusters,
                                                   const reco::PFRecHitDeviceCollection::ConstView inputRecHits,
                                                   const HGCalSoARecHitsExtraDeviceCollection::ConstView inputClusters,
                                                   reco::CaloClusterDeviceCollection::View outputs) const {
    if (numberOfClusters == 0)
      return;

    auto x = make_device_view<float>(queue, outputs.position().x(), numberOfClusters);
    alpaka::memset(queue, x, 0x0);
    auto y = make_device_view<float>(queue, outputs.position().y(), numberOfClusters);
    alpaka::memset(queue, y, 0x0);
    auto z = make_device_view<float>(queue, outputs.position().z(), numberOfClusters);
    alpaka::memset(queue, z, 0x0);
    auto cells = make_device_view<int>(queue, outputs.position().cells(), numberOfClusters);
    alpaka::memset(queue, cells, 0x0);
    auto energy = make_device_view<float>(queue, outputs.energy().energy(), numberOfClusters);
    alpaka::memset(queue, energy, 0x0);

    const auto items = 256u;

    const auto nHits = inputClusters.metadata().size();
    const auto hitsWorkDiv = make_workdiv<Acc1D>(divide_up_by(nHits, items), items);
    alpaka::exec<Acc1D>(queue, hitsWorkDiv, KernelAccumulate<CAL>{}, inputRecHits, inputClusters, outputs);

    const auto clustersWorkDiv = make_workdiv<Acc1D>(divide_up_by(numberOfClusters, items), items);
    alpaka::exec<Acc1D>(queue, clustersWorkDiv, KernelFinalise{}, outputs, numberOfClusters);
  }

  template class BarrelLayerClustersSoAAlgoWrapper<barrelclue::ECAL>;
  template class BarrelLayerClustersSoAAlgoWrapper<barrelclue::HCAL>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
