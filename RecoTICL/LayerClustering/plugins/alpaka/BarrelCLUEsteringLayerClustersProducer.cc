#include <array>
#include <span>

// PFRecHitDeviceCollection pulls in Eigen, which must precede any other SoA header
#include "DataFormats/ParticleFlowReco/interface/alpaka/PFRecHitDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/typedefs.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/SynchronizingEDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

#include "BarrelCLUEsteringAlgoWrapper.h"
#include "BarrelCalorimeters.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  template <typename CAL>
  class BarrelCLUEsteringLayerClustersProducer : public stream::SynchronizingEDProducer<> {
  public:
    BarrelCLUEsteringLayerClustersProducer(edm::ParameterSet const& config)
        : SynchronizingEDProducer(config),
          getTokenDevice_{consumes(config.getParameter<edm::InputTag>("recHits"))},
          getTokenSize_{consumes(config.getParameter<edm::InputTag>("recHits"))},
          noiseToken_{esConsumes(config.getParameter<edm::ESInputTag>("noise"))},
          deviceToken_{produces()},
          deltac_(config.getParameter<double>("deltac")),
          kappa_(config.getParameter<double>("kappa")) {}

    ~BarrelCLUEsteringLayerClustersProducer() override = default;

    void acquire(device::Event const& iEvent, device::EventSetup const& iSetup) override {
      // The PFRecHit SoA is allocated with the capacity of the input CaloRecHits: the actual
      // number of PFRecHits is read from the host product, since the device-side size is not
      // set when there are no input CaloRecHits.
      auto const& deviceInput = iEvent.get(getTokenDevice_);
      nRecHits_ = iEvent.get(getTokenSize_);
      // The per-layer sizes are the CLUEstering batch sizes, which are needed on the host
      algo_.countHitsPerLayer(iEvent.queue(), deviceInput.const_view(), layerSizes_);
    }

    void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {
      auto const& deviceInput = iEvent.get(getTokenDevice_);
      auto const& noise = iSetup.getData(noiseToken_);
      // Allocate output SoA, same size as the input PFRecHits.
      HGCalSoARecHitsExtraDeviceCollection output(iEvent.queue(), nRecHits_);

      algo_.run(iEvent.queue(),
                nRecHits_,
                deltac_,
                kappa_,
                std::span<const uint32_t>(layerSizes_),
                deviceInput.const_view(),
                noise.const_view(),
                output.view());
      iEvent.emplace(deviceToken_, std::move(output));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("recHits", edm::InputTag("TO BE DEFINED"));
      desc.add<edm::ESInputTag>("noise", edm::ESInputTag())
          ->setComment(
              "ECAL: PFRecHitECALParams; HCAL: PFRecHitHCALTopology, which must be produced with "
              "usePFThresholdsFromDB = True");
      desc.add<double>("deltac", 0.0175);
      desc.add<double>("kappa", 1.34);
      // Kept for compatibility with the BarrelLayerClusterProducer configuration
      desc.add<double>("outlierDeltaFactor", 2.)
          ->setComment(
              "Unused: as in the legacy BarrelCLUEAlgo, rechits below the density threshold are only "
              "attached to a higher-density rechit within deltac");
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    device::EDGetToken<reco::PFRecHitDeviceCollection> const getTokenDevice_;
    edm::EDGetTokenT<cms_uint32_t> const getTokenSize_;
    device::ESGetToken<typename CAL::NoiseType, typename CAL::NoiseRecordType> const noiseToken_;
    device::EDPutToken<HGCalSoARecHitsExtraDeviceCollection> const deviceToken_;
    BarrelCLUEsteringAlgoWrapper<CAL> algo_;
    const float deltac_;
    const float kappa_;
    uint32_t nRecHits_;
    std::array<uint32_t, CAL::kNumberOfLayers> layerSizes_;
  };

  using BarrelCLUEsteringLayerClustersProducerECAL = BarrelCLUEsteringLayerClustersProducer<barrelclue::ECAL>;
  using BarrelCLUEsteringLayerClustersProducerHCAL = BarrelCLUEsteringLayerClustersProducer<barrelclue::HCAL>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(BarrelCLUEsteringLayerClustersProducerECAL);
DEFINE_FWK_ALPAKA_MODULE(BarrelCLUEsteringLayerClustersProducerHCAL);
