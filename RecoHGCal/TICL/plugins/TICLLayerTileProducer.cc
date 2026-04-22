// Author: Marco Rovere, marco.rovere@cern.ch
// Date: 05/2019
//
#include <memory>  // unique_ptr

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/TilesHost.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

class TICLLayerTileProducer : public edm::stream::EDProducer<edm::stream::WatchRuns> {
public:
  explicit TICLLayerTileProducer(const edm::ParameterSet &ps);
  ~TICLLayerTileProducer() override {}
  void beginRun(edm::Run const &, edm::EventSetup const &) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_HFNose_token_;
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
        consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_HFNose_clusters"));
    produces<ticl::TICLLayerTilesHFNoseHost>();
  } else {
    clusters_token_ = consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_clusters"));
    produces<ticl::TICLLayerTilesHost>();
    produces<ticl::TICLLayerTilesBarrelHost>("ticlLayerTilesBarrel");
  }
}

void TICLLayerTileProducer::beginRun(edm::Run const &, edm::EventSetup const &es) {
  edm::ESHandle<CaloGeometry> geom = es.getHandle(geometry_token_);
  rhtools_.setGeometry(*geom);
}

void TICLLayerTileProducer::produce(edm::Event &evt, const edm::EventSetup &) {
  using Acc = alpaka_serial_sync::Acc1D;
  // std::unique_ptr<ticl::TICLLayerTilesHFNoseHost> resultHFNose;
  // std::unique_ptr<ticl::TICLLayerTilesHost> result;
  // std::unique_ptr<ticl::TICLLayerTilesBarrelHost> resultBarrel;

  edm::Handle<std::vector<reco::CaloCluster>> cluster_h;
  if (doNose_)
    evt.getByToken(clusters_HFNose_token_, cluster_h);
  else
    evt.getByToken(clusters_token_, cluster_h);

  const auto &layerClusters = *cluster_h;

  std::array<std::vector<float>, ticl::TICLLayerTilesHost::TilesType::nLayers> etas;
  std::array<std::vector<float>, ticl::TICLLayerTilesHost::TilesType::nLayers> phis;
  std::array<std::vector<uint32_t>, ticl::TICLLayerTilesHost::TilesType::nLayers> lcIds;
  std::array<std::vector<float>, ticl::TICLLayerTilesBarrelHost::TilesType::nLayers> barrel_etas;
  std::array<std::vector<float>, ticl::TICLLayerTilesBarrelHost::TilesType::nLayers> barrel_phis;
  std::array<std::vector<uint32_t>, ticl::TICLLayerTilesBarrelHost::TilesType::nLayers> barrel_lcIds;
  // std::array<std::vector<float>, ticl::TICLLayerTilesHFNoseHost::TilesType::nLayers> nose_etas;
  // std::array<std::vector<float>, ticl::TICLLayerTilesHFNoseHost::TilesType::nLayers> nose_phis;
  // std::array<std::vector<uint32_t>, ticl::TICLLayerTilesHFNoseHost::TilesType::nLayers> nose_lcIds;
  auto lcId = 0;
  for (auto const &lc : layerClusters) {
    const auto firstHitDetId = lc.hitsAndFractions()[0].first;
    const auto layer = rhtools_.getLayerWithOffset(firstHitDetId);
    bool isBarrelLC = rhtools_.isBarrel(firstHitDetId);
    // if (!isBarrelLC) {
    //   layer += rhtools_.lastLayer(doNose_) * ((rhtools_.zside(firstHitDetId) + 1) >> 1) - 1;
    // }
    // assert(layer >= 0);

    if (isBarrelLC) {
      barrel_etas[layer].push_back(lc.eta());
      barrel_phis[layer].push_back(lc.phi());
      barrel_lcIds[layer].push_back(lcId);
    } else {
      etas[layer].push_back(lc.eta());
      phis[layer].push_back(lc.phi());
      lcIds[layer].push_back(lcId);
    }
    ++lcId;
  }

  alpaka_serial_sync::Queue queue(cms::alpakatools::host());
  if (doNose_) {
    std::array<int, ticl::TICLLayerTilesHFNoseHost::TilesType::nLayers> nose_sizes;
    std::transform(etas.begin(),
                   etas.begin() + ticl::TICLLayerTilesHFNoseHost::TilesType::nLayers,
                   nose_sizes.begin(),
                   [](const auto &layer) { return layer.size(); });

    auto resultHFNose = std::make_unique<ticl::TICLLayerTilesHFNoseHost>(nose_sizes);
    for (auto layer = 0; layer < ticl::TICLLayerTilesHFNoseHost::TilesType::nLayers; ++layer) {
      (*resultHFNose)[layer].template fill<Acc>(queue, etas[layer], phis[layer], lcIds[layer]);
    }
    evt.put(std::move(resultHFNose));
  } else {
    std::array<int, ticl::TICLLayerTilesBarrelHost::TilesType::nLayers> barrel_sizes;
    std::array<int, ticl::TICLLayerTilesHost::TilesType::nLayers> sizes;
    std::ranges::transform(etas, sizes.begin(), [](const auto &layer) { return layer.size(); });
    std::ranges::transform(barrel_etas, barrel_sizes.begin(), [](const auto &layer) { return layer.size(); });

    auto resultBarrel = std::make_unique<ticl::TICLLayerTilesBarrelHost>(barrel_sizes);
    auto result = std::make_unique<ticl::TICLLayerTilesHost>(sizes);

    for (auto layer = 0; layer < ticl::TICLLayerTilesHost::TilesType::nLayers; ++layer) {
      (*result)[layer].template fill<Acc>(queue, etas[layer], phis[layer], lcIds[layer]);
    }
    for (auto barrel_layer = 0; barrel_layer < ticl::TICLLayerTilesBarrelHost::TilesType::nLayers; ++barrel_layer) {
      (*resultBarrel)[barrel_layer].template fill<Acc>(
          queue, barrel_etas[barrel_layer], barrel_phis[barrel_layer], barrel_lcIds[barrel_layer]);
    }
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
