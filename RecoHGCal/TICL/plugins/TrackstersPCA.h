#ifndef RECOHGCAL_TICL_TRACKSTERSPCA_H
#define RECOHGCAL_TICL_TRACKSTERSPCA_H

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include <vector>

namespace ticl {
  void assignPCAtoTracksters(std::vector<Trackster> &,
                             const std::vector<reco::CaloCluster> &,
                             const edm::ValueMap<std::pair<float, float>> &,
                             hgcal::RecHitTools rhtools,
                             double,
                             bool energyWeight = true,
                             bool cleaning_selection = false);
}
#endif
