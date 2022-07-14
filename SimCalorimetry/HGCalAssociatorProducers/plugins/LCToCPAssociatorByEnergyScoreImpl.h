// Original Author: Marco Rovere

#include <vector>
#include <map>
#include <unordered_map>
#include <memory>  // shared_ptr

#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "SimDataFormats/Associations/interface/LayerClusterToCaloParticleAssociator.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "SimCalorimetry/HGCalAssociatorProducers/interface/AssociatorTools.h"

namespace edm {
  class EDProductGetter;
}

namespace hgcal {
  // This structure is used both for LayerClusters and CaloParticles storing their id and the fraction of a hit
  // that belongs to the LayerCluster or CaloParticle. The meaning of the operator is extremely important since
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

  // This object connects a LayerCluster, identified through its id (lcId), with a vector of pairs containing all the CaloParticles
  // (via their ids (cpIds)) that share at least one cell with the LayerCluster. In that pair it
  // stores the score (lcId->(cpId,score)). Keep in mind that the association is not unique, since there could be several instances
  // of the same CaloParticle from several related SimClusters that each contributed to the same LayerCluster.
  typedef std::vector<std::vector<std::pair<unsigned int, float>>> layerClusterToCaloParticle;
  // This is used to save the caloParticleOnLayer structure for all CaloParticles in each layer.
  // It is not exactly what is returned outside, but out of its entries, the output object is build.
  typedef std::vector<std::vector<hgcal::caloParticleOnLayer>> caloParticleToLayerCluster;
  //This is the output of the makeConnections function that contain all the work with CP2LC and LC2CP
  //association. It will be read by the relevant associateSimToReco and associateRecoToSim functions to
  //provide the final product.
  typedef std::tuple<layerClusterToCaloParticle, caloParticleToLayerCluster> association;
}  // namespace hgcal

class LCToCPAssociatorByEnergyScoreImpl : public hgcal::LayerClusterToCaloParticleAssociatorBaseImpl {
public:
  explicit LCToCPAssociatorByEnergyScoreImpl(edm::EDProductGetter const &,
                                             bool,
                                             std::shared_ptr<hgcal::RecHitTools>,
                                             const std::unordered_map<DetId, const HGCRecHit *> *);

  hgcal::RecoToSimCollection associateRecoToSim(const edm::Handle<reco::CaloClusterCollection> &cCH,
                                                const edm::Handle<CaloParticleCollection> &cPCH) const override;

  hgcal::SimToRecoCollection associateSimToReco(const edm::Handle<reco::CaloClusterCollection> &cCH,
                                                const edm::Handle<CaloParticleCollection> &cPCH) const override;

private:
  const bool hardScatterOnly_;
  std::shared_ptr<hgcal::RecHitTools> recHitTools_;
  const std::unordered_map<DetId, const HGCRecHit *> *hitMap_;
  unsigned layers_;
  edm::EDProductGetter const *productGetter_;
  hgcal::association makeConnections(const edm::Handle<reco::CaloClusterCollection> &cCH,
                                     const edm::Handle<CaloParticleCollection> &cPCH) const;
};
