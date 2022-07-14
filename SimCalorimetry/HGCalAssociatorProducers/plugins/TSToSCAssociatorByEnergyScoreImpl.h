// Original Author: Leonardo Cristella

#include <vector>
#include <map>
#include <unordered_map>
#include <memory>  // shared_ptr

#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "SimDataFormats/Associations/interface/TracksterToSimClusterAssociator.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "SimCalorimetry/HGCalAssociatorProducers/interface/AssociatorTools.h"

namespace edm {
  class EDProductGetter;
}

namespace hgcal {
  struct detIdInfoInCluster {
    bool operator==(const detIdInfoInCluster &o) const { return clusterId == o.clusterId; };
    long unsigned int clusterId;
    float fraction;
    detIdInfoInCluster(long unsigned int cId, float fr) {
      clusterId = cId;
      fraction = fr;
    }
  };

  typedef std::vector<std::vector<std::pair<unsigned int, float>>> tracksterToSimCluster;
  typedef std::vector<hgcal::simClusterOnLayer> simClusterToTrackster;
  typedef std::tuple<tracksterToSimCluster, simClusterToTrackster> association;
}  // namespace hgcal

class TSToSCAssociatorByEnergyScoreImpl : public hgcal::TracksterToSimClusterAssociatorBaseImpl {
public:
  explicit TSToSCAssociatorByEnergyScoreImpl(edm::EDProductGetter const &,
                                             bool,
                                             std::shared_ptr<hgcal::RecHitTools>,
                                             const std::unordered_map<DetId, const HGCRecHit *> *);

  hgcal::RecoToSimCollectionTracksters associateRecoToSim(const edm::Handle<ticl::TracksterCollection> &tCH,
                                                          const edm::Handle<reco::CaloClusterCollection> &lCCH,
                                                          const edm::Handle<SimClusterCollection> &sCCH) const override;

  hgcal::SimToRecoCollectionTracksters associateSimToReco(const edm::Handle<ticl::TracksterCollection> &tCH,
                                                          const edm::Handle<reco::CaloClusterCollection> &lCCH,
                                                          const edm::Handle<SimClusterCollection> &sCCH) const override;

private:
  const bool hardScatterOnly_;
  std::shared_ptr<hgcal::RecHitTools> recHitTools_;
  const std::unordered_map<DetId, const HGCRecHit *> *hitMap_;
  unsigned layers_;
  edm::EDProductGetter const *productGetter_;
  hgcal::association makeConnections(const edm::Handle<ticl::TracksterCollection> &tCH,
                                     const edm::Handle<reco::CaloClusterCollection> &lCCH,
                                     const edm::Handle<SimClusterCollection> &sCCH) const;
};
