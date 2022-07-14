// Original Author: Leonardo Cristella

#include <vector>
#include <map>
#include <unordered_map>
#include <memory>  // shared_ptr

#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "SimDataFormats/Associations/interface/LayerClusterToSimClusterAssociator.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "SimCalorimetry/HGCalAssociatorProducers/interface/AssociatorTools.h"

namespace edm {
  class EDProductGetter;
}

namespace hgcal {
  // This structure is used both for LayerClusters and SimClusters storing their id and the fraction of a hit
  // that belongs to the LayerCluster or SimCluster. The meaning of the operator is extremely important since
  // this struct will be used inside maps and other containers and when searching for one particular occurence
  // only the clusterId member will be used in the check, skipping the fraction part.
  struct detIdInfoInCluster {
    bool operator==(const detIdInfoInCluster &o) const { return clusterId == o.clusterId; };
    long unsigned int clusterId;
    float fraction;
    detIdInfoInCluster(long unsigned int cId, float fr) {
      clusterId = cId;
      fraction = fr;
    }
  };

  // This object connects a LayerCluster, identified through its id (lcId), with a vector of pairs containing all the SimClusters
  // (via their ids (scIds)) that share at least one cell with the LayerCluster. In that pair it
  // stores the score (lcId->(scId,score)).
  typedef std::vector<std::vector<std::pair<unsigned int, float>>> layerClusterToSimCluster;
  // This is used to save the simClusterOnLayer structure for all simClusters in each layer.
  // It is not exactly what is returned outside, but out of its entries, the output object is build.
  typedef std::vector<std::vector<hgcal::simClusterOnLayer>> simClusterToLayerCluster;
  //This is the output of the makeConnections function that contain all the work with SC2LC and LC2SC
  //association. It will be read by the relevant associateSimToReco and associateRecoToSim functions to
  //provide the final product.
  typedef std::tuple<layerClusterToSimCluster, simClusterToLayerCluster> association;
}  // namespace hgcal

class LCToSCAssociatorByEnergyScoreImpl : public hgcal::LayerClusterToSimClusterAssociatorBaseImpl {
public:
  explicit LCToSCAssociatorByEnergyScoreImpl(edm::EDProductGetter const &,
                                             bool,
                                             std::shared_ptr<hgcal::RecHitTools>,
                                             const std::unordered_map<DetId, const HGCRecHit *> *);

  hgcal::RecoToSimCollectionWithSimClusters associateRecoToSim(
      const edm::Handle<reco::CaloClusterCollection> &cCH,
      const edm::Handle<SimClusterCollection> &sCCH) const override;

  hgcal::SimToRecoCollectionWithSimClusters associateSimToReco(
      const edm::Handle<reco::CaloClusterCollection> &cCH,
      const edm::Handle<SimClusterCollection> &sCCH) const override;

private:
  const bool hardScatterOnly_;
  std::shared_ptr<hgcal::RecHitTools> recHitTools_;
  const std::unordered_map<DetId, const HGCRecHit *> *hitMap_;
  unsigned layers_;
  edm::EDProductGetter const *productGetter_;
  hgcal::association makeConnections(const edm::Handle<reco::CaloClusterCollection> &cCH,
                                     const edm::Handle<SimClusterCollection> &sCCH) const;
};
