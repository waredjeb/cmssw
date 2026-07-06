#include <alpaka/alpaka.hpp>

#include "DataFormats/HGCalReco/interface/alpaka/TracksterSoADeviceCollection.h"
#include "DataFormats/HGCalReco/interface/TICLGraph.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "PhysicsTools/PyTorchAlpaka/interface/TensorCollection.h"
#include "PhysicsTools/PyTorchAlpaka/interface/alpaka/AlpakaModel.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class TracksterLinkingByGNNProducer : public stream::EDProducer<> {
  public:
    TracksterLinkingByGNNProducer(const edm::ParameterSet &params);

    void produce(device::Event &event, const device::EventSetup &event_setup) override;
    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  private:
    const device::EDGetToken<TrackstersSoADeviceCollection> inputs_token_;
    const device::EDPutToken<TrackstersGNNOutputSoADeviceCollection> outputs_token_;
	torch::AlpakaModel model_;
  };

  TracksterLinkingByGNNProducer::TracksterLinkingByGNNProducer(edm::ParameterSet const &params)
      : EDProducer<>(params),
        inputs_token_{consumes(params.getParameter<edm::InputTag>("inputs"))},
		outputs_token_{produces()},
		model_(params.getParameter<edm::FileInPath>("model").fullPath()) {}

  void TracksterLinkingByGNNProducer::produce(device::Event &event, const device::EventSetup &event_setup) {

    // get data
    auto &inputs = const_cast<TrackstersSoADeviceCollection &>(event.get(inputs_token_));
    const size_t numNodes = inputs.const_view<GNNNodeSoA>().metadata().size();
    const size_t numEdges = inputs.const_view<GNNEdgeSoA>().metadata().size();
    auto outputs = TrackstersGNNOutputSoADeviceCollection(numEdges, event.queue());
    outputs.zeroInitialise(event.queue());

    if (numNodes > 0) {
      // metadata for automatic tensor conversion
      auto node_records = inputs.view<GNNNodeSoA>().records();
      auto edge_feature_records = inputs.view<GNNEdgeSoA>().records();
      auto edge_index_records = inputs.view<GNNEdgeIndexSoA>().records();
      auto output_records = outputs.view().records();
	  cms::torch::alpakatools::TensorCollection<Queue> inputs_collection(numNodes);

      // Converter can also do full SoA
	  inputs_collection.add<GNNNodeSoA>("nodes",
                                   node_records.barycenter_x(),
                                   node_records.barycenter_y(),
                                   node_records.barycenter_z(),
                                   node_records.barycenter_eta(),
                                   node_records.barycenter_phi(),
                                   node_records.eigenvector0_x(),
                                   node_records.eigenvector0_y(),
                                   node_records.eigenvector0_z(),
                                   node_records.eigenvalue1(),
                                   node_records.eigenvalue2(),
                                   node_records.eigenvalue3(),
                                   node_records.sigmasPCA1(),
                                   node_records.sigmasPCA2(),
                                   node_records.sigmasPCA3(),
                                   node_records.num_LCs(),
                                   node_records.num_hits(),
                                   node_records.raw_energy(),
                                   node_records.raw_em_energy(),
                                   node_records.photon_prob(),
                                   node_records.electron_prob(),
                                   node_records.muon_prob(),
                                   node_records.neutral_pion_prob(),
                                   node_records.charged_hadron_prob(),
                                   node_records.neutral_hadron_prob(),
                                   node_records.z_min(),
                                   node_records.z_max(),
                                   node_records.LC_density(),
                                   node_records.trackster_density(),
                                   node_records.time());

	  inputs_collection.add<GNNEdgeSoA>("edge_features",
                                               numEdges,
                                               edge_feature_records.raw_energy(),
                                               edge_feature_records.barycenter_z(),
                                               edge_feature_records.barycenter_xy(),
                                               edge_feature_records.eigenvector0(),
                                               edge_feature_records.time());

	  inputs_collection.add<GNNEdgeIndexSoA>(
          "edge_index", numEdges, edge_index_records.in(), edge_index_records.out());

	  cms::torch::alpakatools::TensorCollection<Queue> outputs_collection(numEdges);
      outputs_collection.add<GNNOutputSoA>("preds", output_records.score());
      model_.forward(event.queue(), inputs_collection, outputs_collection);

    }
    event.emplace(outputs_token_, std::move(outputs));
    alpaka::wait(event.queue());
  }

  void TracksterLinkingByGNNProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("inputs");
    desc.add<edm::FileInPath>("modelPath");
    descriptions.addWithDefaultLabel(desc);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(TracksterLinkingByGNNProducer);
