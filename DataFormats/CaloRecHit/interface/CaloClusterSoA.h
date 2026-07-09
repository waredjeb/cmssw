
#ifndef DataFormats_CaloRecHit_CaloClusterSoA_H
#define DataFormats_CaloRecHit_CaloClusterSoA_H

// Authors: Simone Balducci, Felice Pantaleo, Wahid Redjeb, Aurora Perego, Leonardo Beltrame

#include "DataFormats/SoATemplate/interface/SoABlocks.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/DetId/interface/DetId.h"
#include <xtd/xtd.h>

namespace reco {

  // clang-format off
  GENERATE_SOA_LAYOUT(CaloClusterSoAPosition,
                      SOA_COLUMN(float, x),
                      SOA_COLUMN(float, y),
                      SOA_COLUMN(float, z),
                      SOA_COLUMN(int, layer),
                      SOA_COLUMN(int, cells),
                      SOA_CONST_ELEMENT_METHODS(
                        SOA_HOST_DEVICE auto r() const { return xtd::sqrt(x() * x() + y() * y()); }
                        SOA_HOST_DEVICE auto phi() const { return xtd::atan2(y(), x()); }
                        SOA_HOST_DEVICE auto eta() const { return xtd::asinh(z() / r()); }
                      ))
  // clang-format on

  GENERATE_SOA_LAYOUT(CaloClusterSoAEnergy,
                      SOA_COLUMN(float, energy),
                      SOA_COLUMN(float, correctedEnergy),
                      SOA_COLUMN(float, correctedEnergyUncertainty))

  GENERATE_SOA_LAYOUT(CaloClusterSoAIndexes,
                      SOA_COLUMN(CaloID, caloID),
                      SOA_COLUMN(CaloCluster::AlgoID, algoID),
                      SOA_COLUMN(DetId, seedID),
                      SOA_COLUMN(uint32_t, flags))

  GENERATE_SOA_LAYOUT(CaloClusterSoATiming, SOA_COLUMN(float, time), SOA_COLUMN(float, timeError))

  GENERATE_SOA_BLOCKS(CaloClusterSoALayout,
                      SOA_BLOCK(position, CaloClusterSoAPosition),
                      SOA_BLOCK(energy, CaloClusterSoAEnergy),
                      SOA_BLOCK(indexes, CaloClusterSoAIndexes),
                      SOA_BLOCK(timing, CaloClusterSoATiming))

  using CaloClusterSoA = CaloClusterSoALayout<>;
  using CaloClusterSoAView = CaloClusterSoA::View;
  using CaloClusterSoAConstView = CaloClusterSoA::ConstView;

}  // namespace reco

#endif
