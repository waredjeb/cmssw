// Author: Marco Rovere, marco.rovere@cern.ch
// Date: 05/2019
//
#include <cassert>

#include <memory>  // unique_ptr

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterHostCollection.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

class TICLLayerTileProducer : public edm::stream::EDProducer<edm::stream::WatchRuns> {
public:
  explicit TICLLayerTileProducer(const edm::ParameterSet &ps);
  ~TICLLayerTileProducer() override {}
  void beginRun(edm::Run const &, edm::EventSetup const &) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  edm::EDGetTokenT<reco::CaloClusterHostCollection> clusters_token_;
  edm::EDGetTokenT<reco::CaloClusterHostCollection> clusters_HFNose_token_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometry_token_;
  hgcal::RecHitTools rhtools_;
  std::string detector_;
  bool doNose_;
  bool doBarrel_;
};

TICLLayerTileProducer::TICLLayerTileProducer(const edm::ParameterSet &ps)
    : detector_(ps.getParameter<std::string>("detector")) {
  geometry_token_ = esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>();

  doNose_ = (detector_ == "HFNose");
  doBarrel_ = (detector_ == "Barrel");

  if (doNose_) {
    clusters_HFNose_token_ = consumes(ps.getParameter<edm::InputTag>("layer_HFNose_clusters"));
    produces<TICLLayerTilesHFNose>();
  } else {
    if (doBarrel_) {
      produces<TICLLayerTilesBarrel>("ticlLayerTilesBarrel");
    }
    clusters_token_ = consumes(ps.getParameter<edm::InputTag>("layer_clusters"));
    produces<TICLLayerTiles>();
  }
}

void TICLLayerTileProducer::beginRun(edm::Run const &, edm::EventSetup const &es) {
  edm::ESHandle<CaloGeometry> geom = es.getHandle(geometry_token_);
  rhtools_.setGeometry(*geom);
}

void TICLLayerTileProducer::produce(edm::Event &evt, const edm::EventSetup &) {
  std::unique_ptr<TICLLayerTilesHFNose> resultHFNose;
  std::unique_ptr<TICLLayerTiles> result;
  std::unique_ptr<TICLLayerTilesBarrel> resultBarrel;
  if (doNose_) {
    resultHFNose = std::make_unique<TICLLayerTilesHFNose>();
  } else {
    if (doBarrel_)
      resultBarrel = std::make_unique<TICLLayerTilesBarrel>();
    result = std::make_unique<TICLLayerTiles>();
  }

  const auto &layerClusters = evt.get(doNose_ ? clusters_HFNose_token_ : clusters_token_).const_view();

  const int numberOfClusters = layerClusters.position().metadata().size();
  for (int lcId = 0; lcId < numberOfClusters; ++lcId) {
    const DetId seedId = layerClusters.indexes()[lcId].seedID();
    int layer = rhtools_.getLayerWithOffset(seedId);
    bool isBarrelLC = rhtools_.isBarrel(seedId);
    if (!isBarrelLC) {
      layer += rhtools_.lastLayer(doNose_) * ((rhtools_.zside(seedId) + 1) >> 1) - 1;
    }
    assert(layer >= 0);

    const auto eta = layerClusters.eta(lcId);
    const auto phi = layerClusters.phi(lcId);

    if (doNose_) {
      resultHFNose->fill(layer, eta, phi, lcId);
      LogDebug("TICLLayerTileProducer") << "Adding layerClusterId: " << lcId << " into bin [eta,phi]: [ "
                                        << (*resultHFNose)[layer].etaBin(eta) << ", "
                                        << (*resultHFNose)[layer].phiBin(phi) << "] for layer: " << layer;
    } else if (doBarrel_ && isBarrelLC) {
      resultBarrel->fill(layer, eta, phi, lcId);
      LogDebug("TICLLayerTileProducer") << "Adding layerClusterId: " << lcId << " into bin [eta,phi]: [ "
                                        << (*resultBarrel)[layer].etaBin(eta) << ", "
                                        << (*resultBarrel)[layer].phiBin(phi) << "] for layer: " << layer;
    } else if (!isBarrelLC) {
      result->fill(layer, eta, phi, lcId);
      LogDebug("TICLLayerTileProducer") << "Adding layerClusterId: " << lcId << " into bin [eta,phi]: [ "
                                        << (*result)[layer].etaBin(eta) << ", "
                                        << (*result)[layer].phiBin(phi) << "] for layer: " << layer;
    }
  }

  if (doNose_)
    evt.put(std::move(resultHFNose));
  else {
    if (doBarrel_)
      evt.put(std::move(resultBarrel), "ticlLayerTilesBarrel");
    else
      evt.put(std::move(result));
  }
}

void TICLLayerTileProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("detector", "HGCAL");
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>("layer_HFNose_clusters", edm::InputTag("hgcalLayerClustersHFNose"));
  descriptions.add("ticlLayerTileProducer", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TICLLayerTileProducer);
