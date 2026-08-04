// Benchmark-only: extract the few HGCalSoARecHits columns that the layer
// clustering actually reads (dim1, dim2, layer, energy, sigmaNoise, detid) into
// plain std::vector products, plus forward layerSizes. Plain vectors deserialize
// far faster than the packed 12-column SoA, so a clustering-only benchmark that
// replays these is limited by clustering, not by reading the rechits.
//
// It copies the device SoA to a host mirror first, so it works on any backend
// (staging can run on CUDA). Not part of the reconstruction chain.
#include <cstdint>
#include <vector>

#include <alpaka/alpaka.hpp>

#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class HGCalSoARecHitsSlimExtractor : public stream::EDProducer<> {
  public:
    HGCalSoARecHitsSlimExtractor(edm::ParameterSet const& config)
        : EDProducer(config),
          soaToken_(consumes(config.getParameter<edm::InputTag>("src"))),
          layerSizesToken_(consumes<std::vector<uint32_t>>(
              edm::InputTag(config.getParameter<edm::InputTag>("src").label(), "layerSizes"))),
          dim1Token_{produces("dim1")},
          dim2Token_{produces("dim2")},
          energyToken_{produces("energy")},
          sigmaToken_{produces("sigmaNoise")},
          layerToken_{produces("layer")},
          detidToken_{produces("detid")},
          layerSizesOut_{produces("layerSizes")} {}

    void produce(device::Event& iEvent, device::EventSetup const&) override {
      auto const& soa = iEvent.get(soaToken_);
      const int n = soa.const_view().metadata().size();
      // Copy to a host mirror so the columns can be read on the host regardless
      // of backend (a no-op copy on the serial backend).
      HGCalSoARecHitsHostCollection hostSoA{iEvent.queue(), n};
      alpaka::memcpy(iEvent.queue(), hostSoA.buffer(), soa.const_buffer());
      alpaka::wait(iEvent.queue());
      auto const v = hostSoA.const_view();
      std::vector<float> dim1(n), dim2(n), energy(n), sigma(n);
      std::vector<int32_t> layer(n);
      std::vector<uint32_t> detid(n);
      for (int i = 0; i < n; ++i) {
        dim1[i] = v[i].dim1();
        dim2[i] = v[i].dim2();
        energy[i] = v[i].energy();
        sigma[i] = v[i].sigmaNoise();
        layer[i] = v[i].layer();
        detid[i] = v[i].detid();
      }
      iEvent.emplace(dim1Token_, std::move(dim1));
      iEvent.emplace(dim2Token_, std::move(dim2));
      iEvent.emplace(energyToken_, std::move(energy));
      iEvent.emplace(sigmaToken_, std::move(sigma));
      iEvent.emplace(layerToken_, std::move(layer));
      iEvent.emplace(detidToken_, std::move(detid));
      iEvent.emplace(layerSizesOut_, iEvent.get(layerSizesToken_));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("src", edm::InputTag("hltHgcalSoARecHitsProducer"));
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    const device::EDGetToken<HGCalSoARecHitsDeviceCollection> soaToken_;
    const edm::EDGetTokenT<std::vector<uint32_t>> layerSizesToken_;
    const edm::EDPutTokenT<std::vector<float>> dim1Token_;
    const edm::EDPutTokenT<std::vector<float>> dim2Token_;
    const edm::EDPutTokenT<std::vector<float>> energyToken_;
    const edm::EDPutTokenT<std::vector<float>> sigmaToken_;
    const edm::EDPutTokenT<std::vector<int32_t>> layerToken_;
    const edm::EDPutTokenT<std::vector<uint32_t>> detidToken_;
    const edm::EDPutTokenT<std::vector<uint32_t>> layerSizesOut_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(HGCalSoARecHitsSlimExtractor);
