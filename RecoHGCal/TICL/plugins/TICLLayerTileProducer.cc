// Author: Marco Rovere, marco.rovere@cern.ch
// Date: 05/2019
//
#include <memory>  // unique_ptr

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterHostCollection.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

class TICLLayerTileProducer : public edm::stream::EDProducer<> {
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
};

TICLLayerTileProducer::TICLLayerTileProducer(const edm::ParameterSet &ps)
    : detector_(ps.getParameter<std::string>("detector")) {
  geometry_token_ = esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>();

  doNose_ = (detector_ == "HFNose");

  if (doNose_) {
    clusters_HFNose_token_ =
        consumes<reco::CaloClusterHostCollection>(ps.getParameter<edm::InputTag>("layer_HFNose_clusters"));
    produces<TICLLayerTilesHFNose>();
  } else {
    clusters_token_ = consumes<reco::CaloClusterHostCollection>(ps.getParameter<edm::InputTag>("layer_clusters"));
    produces<TICLLayerTiles>();
    produces<TICLLayerTilesBarrel>("ticlLayerTilesBarrel");
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
    resultBarrel = std::make_unique<TICLLayerTilesBarrel>();
    result = std::make_unique<TICLLayerTiles>();
  }

  edm::Handle<reco::CaloClusterHostCollection> cluster_h;
  if (doNose_)
    evt.getByToken(clusters_HFNose_token_, cluster_h);
  else
    evt.getByToken(clusters_token_, cluster_h);

  const auto &layerClusters = *cluster_h;
  int lcId = 0;
  for (auto lc_idx = 0; lc_idx < layerClusters.view().position().metadata().size(); ++lc_idx) {
    auto layer = layerClusters.view().position()[lc_idx].layer();
    assert(layer >= 0);

    const auto seed_detid = layerClusters.view().indexes()[lc_idx].seedID();

    const auto isBarrelLC = rhtools_.isBarrel(seed_detid);
    if (!isBarrelLC) {
      layer += rhtools_.lastLayer(doNose_) * ((rhtools_.zside(seed_detid) + 1) >> 1) - 1;
    }

    if (doNose_) {
      resultHFNose->fill(layer, layerClusters.view().eta(lc_idx), layerClusters.view().phi(lc_idx), lc_idx);
    } else if (isBarrelLC) {
      resultBarrel->fill(layer, layerClusters.view().eta(lc_idx), layerClusters.view().phi(lc_idx), lc_idx);
    } else {
      result->fill(layer, layerClusters.view().eta(lc_idx), layerClusters.view().phi(lc_idx), lc_idx);
    }
    LogDebug("TICLLayerTileProducer") << "Adding layerClusterId: " << lc_idx << " into bin [eta,phi]: [ "
                                      << (*result)[layer].etaBin(layerClusters.view().eta(lc_idx)) << ", "
                                      << (*result)[layer].phiBin(layerClusters.view().phi(lc_idx))
                                      << "] for layer: " << layer << std::endl;
  }
  if (doNose_)
    evt.put(std::move(resultHFNose));
  else {
    evt.put(std::move(resultBarrel), "ticlLayerTilesBarrel");
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
