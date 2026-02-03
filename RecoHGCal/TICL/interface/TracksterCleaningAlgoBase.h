#ifndef RecoHGCal_TICL_TracksterCleaningAlgoBase_H__
#define RecoHGCal_TICL_TracksterCleaningAlgoBase_H__

#include <vector>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

namespace ticl {

class TracksterCleaningAlgoBase {
public:
  TracksterCleaningAlgoBase(const edm::ParameterSet& conf, edm::ConsumesCollector)
      : algo_verbosity_(conf.getParameter<int>("algo_verbosity")) {}
  virtual ~TracksterCleaningAlgoBase() = default;

  struct Inputs {
    const edm::Event& ev;
    const edm::EventSetup& es;
    const std::vector<ticl::Trackster>& linked;
    const std::vector<ticl::Trackster>& clue3d;
    const std::vector<reco::CaloCluster>& layerClusters;
    const std::vector<std::vector<unsigned int>>& map; // indices of tracksters associated to each linkedTrackster
    Inputs(const edm::Event& eV,
           const edm::EventSetup& eS,
           const std::vector<ticl::Trackster>& l,
           const std::vector<ticl::Trackster>& c3d,
           const std::vector<reco::CaloCluster>& lc,
           const std::vector<std::vector<unsigned int>>& m)
      : ev(eV), es(eS), linked(l), clue3d(c3d), layerClusters(lc), map(m) {}
  };

  virtual void cleanTracksters(const Inputs& input,
                               std::vector<ticl::Trackster>& outTracksters,
                               std::vector<std::vector<unsigned int>>& outMap) const = 0;

  virtual void initialize() {}
  virtual void setEvent(edm::Event&, edm::EventSetup const&) {}
  static void fillPSetDescription(edm::ParameterSetDescription& d) { d.add<int>("algo_verbosity",0); }

protected:
  int algo_verbosity_{0};
};


}  // namespace ticl

#endif
