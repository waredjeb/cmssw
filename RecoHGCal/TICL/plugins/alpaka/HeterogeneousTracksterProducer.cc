#include <alpaka/alpaka.hpp>
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "DataFormats/HGCalReco/interface/HGCalSoAClusters.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/TilesHost.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/stringize.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "RecoHGCal/TICL/interface/alpaka/PatternRecognitionAlgoBase.h"
#include "RecoHGCal/TICL/plugins/PatternRecognitionPluginFactory.h"
#include "CLUEstering/CLUEstering.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <memory>
#include <string>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // Forward declaration of wrapper class
  class PatternRecognitionByCLUEsteringWrapper;

  class HeterogeneousTracksterProducer : public stream::EDProducer<> {
  public:
    HeterogeneousTracksterProducer(edm::ParameterSet const& config)
        : EDProducer(config),
          // detector_(config.getParameter<std::string>("detector")),
          // doNose_(detector_ == "HFNose"),
          deviceTokenSoAClusters_{consumes(config.getParameter<edm::InputTag>("layerClusters"))},
          layer_clusters_tiles_token_(consumes<ticl::TICLLayerTilesHost>(config.getParameter<edm::InputTag>("layer_clusters_tiles"))),
          legacyTrackstersToken_{produces()}
  {
      auto plugin = config.getParameter<std::string>("patternRecognitionBy");
      auto pluginPSet = config.getParameter<edm::ParameterSet>("pluginPatternRecognitionBy" + plugin);

      // Construct backend-specific plugin name: "alpaka_serial_sync::CLUEstering" or "alpaka_cuda_async::CLUEstering"
      std::string backendSpecificName = std::string(EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE)) + "::" + plugin;

      // Use backend-independent factory with backend-specific plugin name
      auto algoWrapper = PatternRecognitionFactoryPortable::get()->create(backendSpecificName, pluginPSet);

      // Cast to wrapper to get backend-specific implementation
      auto* wrapper = dynamic_cast<PatternRecognitionByCLUEsteringWrapper*>(algoWrapper.release());
      if (!wrapper) {
        throw cms::Exception("HeterogeneousTracksterProducer")
            << "Failed to cast pattern recognition algorithm to wrapper type";
      }
      algo_.reset(wrapper->getImpl());
    }
    ~HeterogeneousTracksterProducer() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalSoALayerClusters"));
      // Still needed? I guess I still need to save them here and then copy to algo
      // desc.add<double>("rho_c", 0.6);
      // desc.add<std::vector<double>>("dc", {2., 2., 2});
      // desc.add<std::vector<double>>("dm", {1.8, 1.8, 2});l
      desc.add<std::string>("patternRecognitionBy", "CLUEstering");

      edm::ParameterSetDescription pluginDesc;
      // Use backend-independent factory with backend-specific type name
      std::string backendSpecificType = std::string(EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE)) + "::CLUEstering";
      pluginDesc.addNode(edm::PluginDescription<PatternRecognitionFactoryPortable>("type", backendSpecificType, true));
      desc.add<edm::InputTag>("layer_clusters_tiles", edm::InputTag("ticlLayerTileProducer"));
      desc.add<edm::ParameterSetDescription>("pluginPatternRecognitionByCLUEstering", pluginDesc);

      descriptions.addWithDefaultLabel(desc);
    }

    void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {
      const auto& lc = iEvent.get(deviceTokenSoAClusters_);
      const auto& tiles = iEvent.get(layer_clusters_tiles_token_);

      auto tilesView = tiles.view();

      for(int currentLayer = 1; currentLayer < ticl::TICLLayerTilesHost::TilesType::nLayers; currentLayer++){
        auto const layerTile = tilesView[currentLayer];
        for(int currentTile = 0; currentTile < ticl::TICLLayerTilesHost::TilesType::nBins; currentTile++) {
          if(layerTile.count(currentTile) > 0)
            std::cout << "LayerTile on layer " << currentLayer << " Has " << layerTile.count(currentTile) << " For tile " << currentTile << std::endl;
       }
      }
      
      auto tracksters = std::vector<ticl::Trackster>();
      auto& queue = iEvent.queue();
      algo_->makeTracksters(queue, lc, tracksters);

      iEvent.emplace(legacyTrackstersToken_, std::move(tracksters));
    }

  private:
    // std::string detector_;
    // bool doNose_;
    device::EDGetToken<HGCalSoAClustersDeviceCollection> const deviceTokenSoAClusters_;
    edm::EDGetTokenT<ticl::TICLLayerTilesHost> layer_clusters_tiles_token_;
    edm::EDPutTokenT<std::vector<ticl::Trackster>> const legacyTrackstersToken_;
    std::unique_ptr<PatternRecognitionAlgoBase> algo_;
    std::unique_ptr<PatternRecognitionAlgoBase> myAlgoHFNose_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(HeterogeneousTracksterProducer);
