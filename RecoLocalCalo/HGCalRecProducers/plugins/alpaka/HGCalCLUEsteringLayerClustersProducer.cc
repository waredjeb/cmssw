#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include "HGCalCLUEsteringAlgoWrapper.h"

// Heterogeneous (alpaka) producer that runs the external CLUEstering library on
// device to build HGCal 2D per-layer layer-clusters. It is a drop-in
// replacement for HGCalSoARecHitsLayerClustersProducer: it consumes the same
// HGCalSoARecHitsDeviceCollection and produces the same
// HGCalSoARecHitsExtraDeviceCollection, so the rest of the device chain
// (HGCalSoARecHitsProducer upstream, HGCalSoALayerClustersProducer downstream)
// is reused unchanged.
namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class HGCalCLUEsteringLayerClustersProducer : public stream::EDProducer<> {
  public:
    HGCalCLUEsteringLayerClustersProducer(edm::ParameterSet const& config)
        : EDProducer(config),
          getTokenDevice_{consumes(config.getParameter<edm::InputTag>("hgcalRecHitsSoA"))},
          deviceToken_{produces()},
          deltac_((float)config.getParameter<double>("deltac")),
          kappa_((float)config.getParameter<double>("kappa")),
          outlierDeltaFactor_((float)config.getParameter<double>("outlierDeltaFactor")) {}

    ~HGCalCLUEsteringLayerClustersProducer() override = default;

    void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {
      auto const& deviceInput = iEvent.get(getTokenDevice_);
      auto const input_v = deviceInput.view();
      // Allocate output SoA, same size as the input RecHit SoA.
      HGCalSoARecHitsExtraDeviceCollection output(iEvent.queue(), deviceInput->metadata().size());
      auto output_v = output.view();

      algo_.run(
          iEvent.queue(), deviceInput->metadata().size(), deltac_, kappa_, outlierDeltaFactor_, input_v, output_v);
      iEvent.emplace(deviceToken_, std::move(output));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("hgcalRecHitsSoA", edm::InputTag("TO BE DEFINED"));
      desc.add<double>("deltac", 1.3);
      desc.add<double>("kappa", 9.);
      desc.add<double>("outlierDeltaFactor", 2.);
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    device::EDGetToken<HGCalSoARecHitsDeviceCollection> const getTokenDevice_;
    device::EDPutToken<HGCalSoARecHitsExtraDeviceCollection> const deviceToken_;
    HGCalCLUEsteringAlgoWrapper algo_;
    const float deltac_;
    const float kappa_;
    const float outlierDeltaFactor_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(HGCalCLUEsteringLayerClustersProducer);
