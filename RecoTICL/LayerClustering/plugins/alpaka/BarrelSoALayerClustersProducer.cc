// PFRecHitDeviceCollection pulls in Eigen, which must precede any other SoA header
#include "DataFormats/ParticleFlowReco/interface/alpaka/PFRecHitDeviceCollection.h"
#include "DataFormats/CaloRecHit/interface/alpaka/CaloClusterDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/SynchronizingEDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

#include "BarrelCalorimeters.h"
#include "BarrelLayerClustersSoAAlgoWrapper.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  template <typename CAL>
  class BarrelSoALayerClustersProducer : public stream::SynchronizingEDProducer<> {
  public:
    BarrelSoALayerClustersProducer(edm::ParameterSet const& config)
        : SynchronizingEDProducer(config),
          getTokenDeviceRecHits_{consumes(config.getParameter<edm::InputTag>("recHits"))},
          getTokenDeviceClusters_{consumes(config.getParameter<edm::InputTag>("barrelRecHitsLayerClustersSoA"))},
          deviceTokenSoAClusters_{produces()} {}

    ~BarrelSoALayerClustersProducer() override = default;

    void acquire(device::Event const& iEvent, device::EventSetup const& iSetup) override {
      // The number of clusters is needed on the host to allocate the output
      auto const& deviceInputClusters = iEvent.get(getTokenDeviceClusters_);
      auto device_numclusters = cms::alpakatools::make_device_view<const unsigned int>(
          iEvent.queue(), deviceInputClusters.const_view().numberOfClustersScalar());
      auto host_numclusters = cms::alpakatools::make_host_view<unsigned int>(num_clusters_);
      alpaka::memcpy(iEvent.queue(), host_numclusters, device_numclusters);
    }

    void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {
      auto const& deviceInputRecHits = iEvent.get(getTokenDeviceRecHits_);
      auto const& deviceInputClusters = iEvent.get(getTokenDeviceClusters_);

      reco::CaloClusterDeviceCollection output(
          iEvent.queue(), num_clusters_, num_clusters_, num_clusters_, num_clusters_);

      algo_.run(iEvent.queue(),
                num_clusters_,
                deviceInputRecHits.const_view(),
                deviceInputClusters.const_view(),
                output.view());
      iEvent.emplace(deviceTokenSoAClusters_, std::move(output));
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("recHits", edm::InputTag("TO BE DEFINED"));
      desc.add<edm::InputTag>("barrelRecHitsLayerClustersSoA", edm::InputTag("TO BE DEFINED"));
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    device::EDGetToken<reco::PFRecHitDeviceCollection> const getTokenDeviceRecHits_;
    device::EDGetToken<HGCalSoARecHitsExtraDeviceCollection> const getTokenDeviceClusters_;
    device::EDPutToken<reco::CaloClusterDeviceCollection> const deviceTokenSoAClusters_;
    BarrelLayerClustersSoAAlgoWrapper<CAL> algo_;
    unsigned int num_clusters_;
  };

  using BarrelSoALayerClustersProducerECAL = BarrelSoALayerClustersProducer<barrelclue::ECAL>;
  using BarrelSoALayerClustersProducerHCAL = BarrelSoALayerClustersProducer<barrelclue::HCAL>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(BarrelSoALayerClustersProducerECAL);
DEFINE_FWK_ALPAKA_MODULE(BarrelSoALayerClustersProducerHCAL);
