// Benchmark-only: rebuild a HGCalSoARecHits device collection from the slim
// std::vector columns produced by HGCalSoARecHitsSlimExtractor, and forward
// layerSizes. Only the columns the clustering reads are filled; the rest are
// left unspecified (clustering never touches them). This lets the clustering
// step be benchmarked in isolation while the (lightweight) vectors are cheap to
// read. Not part of the reconstruction chain.
#include <cstdint>
#include <vector>

#include <alpaka/alpaka.hpp>

#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class HGCalSlimToDeviceProducer : public stream::EDProducer<> {
  public:
    HGCalSlimToDeviceProducer(edm::ParameterSet const& config)
        : EDProducer(config),
          dim1Token_(consumes<std::vector<float>>(tag(config, "dim1"))),
          dim2Token_(consumes<std::vector<float>>(tag(config, "dim2"))),
          energyToken_(consumes<std::vector<float>>(tag(config, "energy"))),
          sigmaToken_(consumes<std::vector<float>>(tag(config, "sigmaNoise"))),
          layerToken_(consumes<std::vector<int32_t>>(tag(config, "layer"))),
          detidToken_(consumes<std::vector<uint32_t>>(tag(config, "detid"))),
          layerSizesToken_(consumes<std::vector<uint32_t>>(tag(config, "layerSizes"))),
          deviceToken_{produces()},
          layerSizesOut_{produces("layerSizes")} {}

    void produce(device::Event& iEvent, device::EventSetup const&) override {
      auto const& dim1 = iEvent.get(dim1Token_);
      auto const& dim2 = iEvent.get(dim2Token_);
      auto const& energy = iEvent.get(energyToken_);
      auto const& sigma = iEvent.get(sigmaToken_);
      auto const& layer = iEvent.get(layerToken_);
      auto const& detid = iEvent.get(detidToken_);
      const int n = dim1.size();

      HGCalSoARecHitsHostCollection host{iEvent.queue(), n};
      auto hv = host.view();
      for (int i = 0; i < n; ++i) {
        hv[i].dim1() = dim1[i];
        hv[i].dim2() = dim2[i];
        hv[i].energy() = energy[i];
        hv[i].sigmaNoise() = sigma[i];
        hv[i].layer() = layer[i];
        hv[i].detid() = detid[i];
      }
      HGCalSoARecHitsDeviceCollection deviceProduct{iEvent.queue(), n};
      alpaka::memcpy(iEvent.queue(), deviceProduct.buffer(), host.const_buffer());
      iEvent.emplace(deviceToken_, std::move(deviceProduct));
      iEvent.emplace(layerSizesOut_, iEvent.get(layerSizesToken_));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("src", edm::InputTag("slimEE"));
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    static edm::InputTag tag(edm::ParameterSet const& c, const char* inst) {
      return edm::InputTag(c.getParameter<edm::InputTag>("src").label(), inst);
    }
    const edm::EDGetTokenT<std::vector<float>> dim1Token_;
    const edm::EDGetTokenT<std::vector<float>> dim2Token_;
    const edm::EDGetTokenT<std::vector<float>> energyToken_;
    const edm::EDGetTokenT<std::vector<float>> sigmaToken_;
    const edm::EDGetTokenT<std::vector<int32_t>> layerToken_;
    const edm::EDGetTokenT<std::vector<uint32_t>> detidToken_;
    const edm::EDGetTokenT<std::vector<uint32_t>> layerSizesToken_;
    const device::EDPutToken<HGCalSoARecHitsDeviceCollection> deviceToken_;
    const edm::EDPutTokenT<std::vector<uint32_t>> layerSizesOut_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(HGCalSlimToDeviceProducer);
