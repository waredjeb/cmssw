#ifndef RecoHGCal_TICL_TrackstersLinkingAlgoBase_H__
#define RecoHGCal_TICL_TrackstersLinkingAlgoBase_H__

#include <memory>
#include <vector>
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/TrackReco/interface/Track.h"
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
      const std::vector<Trackster>& tracksters;
      const std::vector<reco::Track>& trks;

      Inputs(const edm::Event& eV,
             const edm::EventSetup& eS,
             const std::vector<Trackster>& ts,
             const std::vector<reco::Track>& tracks)
          : ev(eV), es(eS), tracksters(ts), trks(tracks) {}
    };

    virtual void trackstersInfo(const Inputs& input) = 0;
  };
}  // namespace ticl

#endif