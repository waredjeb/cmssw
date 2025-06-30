#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/alpaka/TracksterSoADeviceCollection.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLGraph.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "PhysicsTools/PyTorch/interface/Nvtx.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  class TracksterSoAProducer : public stream::EDProducer<> {
  public:
    TracksterSoAProducer(const edm::ParameterSet &params);

    void produce(device::Event &event, const device::EventSetup &event_setup) override;
    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  private:
    const edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_token_;
    const edm::EDGetTokenT<TICLGraph> ticl_graph_token_;
    const edm::EDGetTokenT<std::vector<reco::CaloCluster>> layer_clusters_token_;
    const uint32_t batch_size_; /**< Size of the batch to be produced. */
    const device::EDPutToken<TrackstersSoADeviceCollection> tracksterSoA_token_; /**< Token to store output data. */
  };

  TracksterSoAProducer::TracksterSoAProducer(edm::ParameterSet const &params)
      : EDProducer<>(params),
        tracksters_token_(consumes<std::vector<ticl::Trackster>>(params.getParameter<edm::InputTag>("tracksters"))),
        ticl_graph_token_(consumes<TICLGraph>(params.getParameter<edm::InputTag>("ticlGraph"))),
        layer_clusters_token_(consumes<std::vector<reco::CaloCluster>>(params.getParameter<edm::InputTag>("layerClusters"))),
        batch_size_(params.getParameter<uint32_t>("batchSize")), 
        tracksterSoA_token_{produces()} {}

  void TracksterSoAProducer::produce(device::Event &event, const device::EventSetup &event_setup) {
    auto t1 = std::chrono::high_resolution_clock::now();
    auto const& ticlGraph = event.get(ticl_graph_token_);
    auto const& layerClusters = event.get(layer_clusters_token_);

    // debug stream usage in concurrently scheduled modules
    std::stringstream msg_stream;
    msg_stream << "Data::produce [E: " << event.id().event() << "]";
    auto msg = msg_stream.str();
    NvtxScopedRange produce_range(msg.c_str());

    // create dummy data
    auto collection = TrackstersSoADeviceCollection(batch_size_, event.queue());
    collection.zeroInitialise(event.queue());
    event.emplace(tracksterSoA_token_, std::move(collection));
    alpaka::wait(event.queue());
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "(Data) E: " << event.id().event() << " OK - "
              << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " us" << std::endl;
    produce_range.end();
  }

  /**
   * @brief Describes the allowed configuration parameters for this module.
   * @param descriptions Configuration description object to populate.
   */
  void TracksterSoAProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracksters", edm::InputTag("ticlTrackstersCLUE3DHigh"));
    desc.add<edm::InputTag>("ticlGraph", edm::InputTag("ticlGraph"));
    desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalMergeLayerClusters"));
    desc.add<uint32_t>("batchSize");
    descriptions.addWithDefaultLabel(desc);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(TracksterSoAProducer);
