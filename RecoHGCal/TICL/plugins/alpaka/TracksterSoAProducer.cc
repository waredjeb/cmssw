#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/alpaka/TracksterSoADeviceCollection.h"
#include "DataFormats/HGCalReco/interface/TracksterSoAHostCollection.h"
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
    float detector_size = (2*(3 - 1.5) * (2 * 47));
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
    auto const& tracksters = event.get(tracksters_token_);
    auto const& layerClusters = event.get(layer_clusters_token_);

    // debug stream usage in concurrently scheduled modules
    std::stringstream msg_stream;
    msg_stream << "TracksterSoAProducer::produce [E: " << event.id().event() << "]";
    auto msg = msg_stream.str();
    NvtxScopedRange produce_range(msg.c_str());

    // create dummy data
    auto hostCollection = TrackstersSoAHostCollection(batch_size_, event.queue());
    auto deviceCollection = TrackstersSoADeviceCollection(batch_size_, event.queue());
    TrackstersSoAView& view = hostCollection.view();

    view.trackster_density() = batch_size_ / detector_size;
    
    for (size_t i = 0; i < batch_size_; i++) {
        view.time()[i] = tracksters[i].time();
        view.raw_energy()[i] = tracksters[i].raw_energy();
        view.raw_em_energy()[i] = tracksters[i].raw_em_energy();

        view.barycenter_x()[i] = tracksters[i].barycenter().x();
        view.barycenter_y()[i] = tracksters[i].barycenter().y();
        view.barycenter_z()[i] = tracksters[i].barycenter().z();
        view.barycenter_eta()[i] = tracksters[i].barycenter().eta();
        view.barycenter_phi()[i] = tracksters[i].barycenter().phi();

        view.eigenvector0_x()[i] = tracksters[i].eigenvectors(0).x();
        view.eigenvector0_y()[i] = tracksters[i].eigenvectors(0).y();
        view.eigenvector0_z()[i] = tracksters[i].eigenvectors(0).z();

        view.eigenvalue1()[i] = tracksters[i].eigenvalues()[0];
        view.eigenvalue2()[i] = tracksters[i].eigenvalues()[1];
        view.eigenvalue3()[i] = tracksters[i].eigenvalues()[2];

        view.sigmasPCA1()[i] = tracksters[i].sigmasPCA()[0];
        view.sigmasPCA2()[i] = tracksters[i].sigmasPCA()[1];
        view.sigmasPCA3()[i] = tracksters[i].sigmasPCA()[2];

        view.photon_prob()[i] = tracksters[i].id_probability(ticl::Trackster::ParticleType::photon);
        view.electron_prob()[i] = tracksters[i].id_probability(ticl::Trackster::ParticleType::electron);
        view.muon_prob()[i] = tracksters[i].id_probability(ticl::Trackster::ParticleType::muon);
        view.neutral_pion_prob()[i] = tracksters[i].id_probability(ticl::Trackster::ParticleType::neutral_pion);
        view.charged_hadron_prob()[i] = tracksters[i].id_probability(ticl::Trackster::ParticleType::charged_hadron);
        view.neutral_hadron_prob()[i] = tracksters[i].id_probability(ticl::Trackster::ParticleType::neutral_hadron);

        view.num_LCs()[i] = tracksters[i].vertices().size();

        int hits = 0;
        float z_min = std::numeric_limits<float>::max();
        float z_max = std::numeric_limits<float>::min();
        for (auto& vertex : tracksters[i].vertices()) {
          auto& cluster = layerClusters[vertex];
          hits += cluster.size();

          if (cluster.z() > z_max)
            z_max = cluster.z();

          if (cluster.z() < z_min)
            z_min = cluster.z();
        }

        view.z_min()[i] = z_min;
        view.z_max()[i] = z_max;
        view.LC_density()[i] = view.num_LCs()[i] / detector_size;
    }

    alpaka::memcpy(event.queue(), deviceCollection.buffer(), hostCollection.buffer());
    alpaka::wait(event.queue());

    event.emplace(tracksterSoA_token_, std::move(deviceCollection));
    alpaka::wait(event.queue());
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "(TracksterSoAProducer) E: " << event.id().event() << " OK - "
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
