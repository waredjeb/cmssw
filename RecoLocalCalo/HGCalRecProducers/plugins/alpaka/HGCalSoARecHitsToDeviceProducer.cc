// Copy a HGCalSoARecHits host collection (e.g. one read back from a ROOT file)
// onto the device, so the device clustering step can consume it. This exists so
// the clustering step can be benchmarked in isolation, replaying pre-staged SoA
// rechits without re-running the (expensive) rechit reconstruction. It is not
// part of the normal reconstruction chain.
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

  class HGCalSoARecHitsToDeviceProducer : public stream::EDProducer<> {
  public:
    HGCalSoARecHitsToDeviceProducer(edm::ParameterSet const& config)
        : EDProducer(config),
          hostToken_(consumes<HGCalSoARecHitsHostCollection>(config.getParameter<edm::InputTag>("src"))),
          layerSizesToken_(consumes<std::vector<uint32_t>>(
              edm::InputTag(config.getParameter<edm::InputTag>("src").label(), "layerSizes"))),
          deviceToken_{produces()},
          layerSizesToken2_{produces("layerSizes")} {}

    void produce(device::Event& iEvent, device::EventSetup const&) override {
      auto const& host = iEvent.get(hostToken_);
      const int size = host.view().metadata().size();
      HGCalSoARecHitsDeviceCollection deviceProduct{iEvent.queue(), size};
      alpaka::memcpy(iEvent.queue(), deviceProduct.buffer(), host.const_buffer());
      iEvent.emplace(deviceToken_, std::move(deviceProduct));
      // Forward the per-layer batch sizes so the clustering step (which derives
      // its layerSizes label from hgcalRecHitsSoA) finds them on this module.
      iEvent.emplace(layerSizesToken2_, iEvent.get(layerSizesToken_));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("src", edm::InputTag("hltHgcalSoARecHitsProducer"));
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    const edm::EDGetTokenT<HGCalSoARecHitsHostCollection> hostToken_;
    const edm::EDGetTokenT<std::vector<uint32_t>> layerSizesToken_;
    const device::EDPutToken<HGCalSoARecHitsDeviceCollection> deviceToken_;
    const edm::EDPutTokenT<std::vector<uint32_t>> layerSizesToken2_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(HGCalSoARecHitsToDeviceProducer);
