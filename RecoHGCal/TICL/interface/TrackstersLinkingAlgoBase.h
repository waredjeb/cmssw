#ifndef RecoHGCal_TICL_TrackstersLinkingAlgoBase_H__
#define RecoHGCal_TICL_TrackstersLinkingAlgoBase_H__

#include <memory>
#include <vector>
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "RecoHGCal/TICL/interface/GlobalCache.h"

namespace edm {
  class Event;
  class EventSetup;
}  // namespace edm

namespace ticl {
  class TrackstersLinkingAlgoBaseT {
  public:
    TrackstersLinkingAlgoBaseT(const edm::ParameterSet& conf, const CacheBase* cache){};
    virtual ~TrackstersLinkingAlgoBaseT(){};

    struct Inputs {
      const edm::Event& ev;
      const edm::EventSetup& es;
      const std::vector<reco::CaloCluster>& layerClusters;
      const std::vector<Trackster>& tracksters;

      Inputs(const edm::Event& eV,
             const edm::EventSetup& eS,
             const std::vector<reco::CaloCluster>& lC,
             const std::vector<Trackster>& ts)
          : ev(eV), es(eS), layerClusters(lC), tracksters(ts) {}
    };
  };
}  // namespace ticl

#endif