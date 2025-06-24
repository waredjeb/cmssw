#ifndef RecoHGCal_TICL_TracksterFilterBase_H__
#define RecoHGCal_TICL_TracksterFilterBase_H__

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "FWCore/PluginManager/interface/PluginFactory.h"

#include <memory>
#include <vector>

namespace edm {
  class ParameterSet;
}

namespace ticl {

class TracksterFilterBase {
public:
  explicit TracksterFilterBase(const edm::ParameterSet &) {}
  virtual ~TracksterFilterBase() {}
  virtual void filter(const std::vector<ticl::Trackster> &tracksters,
                      const std::vector<reco::CaloCluster> &layerClusters,
                      std::vector<float> &trackstersMask,
                      hgcal::RecHitTools &rhtools) const = 0;
};
}



typedef edmplugin::PluginFactory<ticl::TracksterFilterBase*(const edm::ParameterSet &)> TracksterFilterFactory;

#endif
// RecoHGCal_TICL_TracksterFilterBase_H__

