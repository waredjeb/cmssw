#include <alpaka/alpaka.hpp>
#include <torch/torch.h>
#include <torch/script.h>

#include "DataFormats/PyTorchTest/interface/alpaka/Collections.h"
#include "DataFormats/HGCalReco/interface/alpaka/TracksterSoADeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/GNNOutputSoADeviceCollection.h"
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
#include "PhysicsTools/PyTorch/interface/AlpakaConfig.h"
#include "PhysicsTools/PyTorch/interface/Model.h"
#include "PhysicsTools/PyTorch/interface/SoAMetadata.h"
#include "PhysicsTools/PyTorch/interface/Nvtx.h"
//#include "PhysicsTools/PyTorch/plugins/alpaka/Kernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using JitModel = cms::torch::alpaka::Model<cms::torch::alpaka::CompilationType::kJustInTime>;

  class TracksterLinkingByGNNProducer : public stream::EDProducer<> {
  public:
    TracksterLinkingByGNNProducer (const edm::ParameterSet &params);

    void produce(device::Event &event, const device::EventSetup &event_setup) override;
    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  private:
    const device::EDGetToken<TrackstersSoADeviceCollection> inputs_token_;
    const device::EDPutToken<TrackstersGNNOutputSoADeviceCollection> outputs_token_;
//    std::unique_ptr<Kernels> kernels_ = nullptr; /**< Kernel utilities for post-inference validation. */
    std::unique_ptr<JitModel> model_;
  };

  TracksterLinkingByGNNProducer::TracksterLinkingByGNNProducer(edm::ParameterSet const &params)
      : EDProducer<>(params),
        inputs_token_{consumes(params.getParameter<edm::InputTag>("inputs"))},
        outputs_token_{produces()}{
    cms::torch::alpaka::set_threading_guard();
    model_ = std::make_unique<JitModel>(params.getParameter<edm::FileInPath>("modelPath").fullPath());
  }

  void TracksterLinkingByGNNProducer::produce(device::Event &event, const device::EventSetup &event_setup) {
    auto t1 = std::chrono::high_resolution_clock::now();

    // debug stream usage in concurrently scheduled modules
    std::stringstream msg_stream;
    msg_stream << "TracksterLinkingGNN::produce [E: " << event.id().event() << "]";
    auto msg = msg_stream.str();
    NvtxScopedRange produce_range(msg.c_str());

    // guard torch internal operations to not conflict with fw execution scheme
    cms::torch::alpaka::Guard<Queue> guard(event.queue());
    // sanity check
    assert(cms::torch::alpaka::queue_hash(event.queue()) == cms::torch::alpaka::current_stream_hash(event.queue()));

    // get data
    auto &inputs = const_cast<TrackstersSoADeviceCollection&>(event.get(inputs_token_));
    const size_t numNodes = inputs.const_view<GNNNodeSoA>().metadata().size();
    const size_t numEdges = inputs.const_view<GNNEdgeSoA>().metadata().size();
    auto outputs = TrackstersGNNOutputSoADeviceCollection(numEdges, event.queue());

    // metadata for automatic tensor conversion
    auto node_records = inputs.view<GNNNodeSoA>().records();
    auto edge_feature_records = inputs.view<GNNEdgeSoA>().records();
    auto edge_index_records = inputs.view<GNNEdgeIndexSoA>().records();
    auto output_records = outputs.view().records();
    cms::torch::alpaka::SoAMetadata<GNNNodeSoA> inputs_metadata(numNodes);

    // Converter can also do full SoA
    inputs_metadata.append_block("nodes", 
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

      inputs_metadata.append_block<GNNEdgeSoA>("edge_features", numEdges,
                            edge_feature_records.raw_energy(),
                            edge_feature_records.barycenter_z(),
                            edge_feature_records.barycenter_xy(),
                            edge_feature_records.eigenvector0(),
                            edge_feature_records.time());

      inputs_metadata.append_block<GNNEdgeIndexSoA>("edge_index", numEdges,
                            edge_index_records.in(),
                            edge_index_records.out());

    cms::torch::alpaka::SoAMetadata<GNNOutputSoA> outputs_metadata(numEdges);
    outputs_metadata.append_block("preds", output_records.score());

    cms::torch::alpaka::ModelMetadata<GNNNodeSoA, GNNOutputSoA> metadata(
        inputs_metadata, outputs_metadata);

    // inference
   NvtxScopedRange move_to_device("Classifier::move_to_device");
   if (cms::torch::alpaka::device(event.queue()) != model_->device()) {
     std::cout << "(TracksterLinkingGNN) E: " << event.id().event() << " Model: " << model_->device() << " -> "
               << cms::torch::alpaka::device(event.queue()) << std::endl;
     model_->to(event.queue());
   }
   assert(cms::torch::alpaka::device(event.queue()) == model_->device());
   move_to_device.end();
    NvtxScopedRange infer_range("TracksterLinkingGNN::inference");
    model_->forward(metadata);
    infer_range.end();

    //kernels_->AssertClassification(event.queue(), outputs);
    alpaka::wait(event.queue());
    event.emplace(outputs_token_, std::move(outputs));
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "(TracksterLinkingGNN) E: " << event.id().event() << " OK - "
              << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << " us" << std::endl;
    produce_range.end();
  }

  void TracksterLinkingByGNNProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("inputs");
    desc.add<edm::FileInPath>("modelPath");
    descriptions.addWithDefaultLabel(desc);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(TracksterLinkingByGNNProducer);
