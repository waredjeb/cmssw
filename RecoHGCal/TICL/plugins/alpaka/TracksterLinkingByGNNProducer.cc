#include <alpaka/alpaka.hpp>
#include <torch/torch.h>
#include <torch/script.h>

#include "DataFormats/PyTorchTest/interface/alpaka/Collections.h"
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
    const device::EDPutToken<torchportable::ClassificationCollection> outputs_token_;
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
    msg_stream << "ClassifierAot::produce [E: " << event.id().event() << "]";
    auto msg = msg_stream.str();
    NvtxScopedRange produce_range(msg.c_str());

    // guard torch internal operations to not conflict with fw execution scheme
    cms::torch::alpaka::Guard<Queue> guard(event.queue());
    // sanity check
    assert(cms::torch::alpaka::queue_hash(event.queue()) == cms::torch::alpaka::current_stream_hash(event.queue()));

    // get data
    auto &inputs = const_cast<TrackstersSoADeviceCollection&>(event.get(inputs_token_));
    const size_t batch_size = inputs.const_view().metadata().size();
    auto outputs = torchportable::ClassificationCollection(batch_size, event.queue());

    // metadata for automatic tensor conversion
    auto input_records = inputs.view().records();
    auto output_records = outputs.view().records();
    cms::torch::alpaka::SoAMetadata<TrackstersSoA> inputs_metadata(batch_size);
    inputs_metadata.append_block("features", input_records.barycenter_x(), input_records.barycenter_y(), input_records.barycenter_z());
    cms::torch::alpaka::SoAMetadata<torchportable::ClassificationSoA> outputs_metadata(batch_size);
    outputs_metadata.append_block("preds", output_records.c1(), output_records.c2());
    cms::torch::alpaka::ModelMetadata<TrackstersSoA, torchportable::ClassificationSoA> metadata(
        inputs_metadata, outputs_metadata);

    // inference
    NvtxScopedRange move_to_device("Classifier::move_to_device");
    if (cms::torch::alpaka::device(event.queue()) != model_->device()) {
      std::cout << "(ClassifierAot) E: " << event.id().event() << " Model: " << model_->device() << " -> "
                << cms::torch::alpaka::device(event.queue()) << std::endl;
      model_->to(event.queue());
    }
    assert(cms::torch::alpaka::device(event.queue()) == model_->device());
    move_to_device.end();
    NvtxScopedRange infer_range("Classifier::inference");
    model_->forward(metadata);
    infer_range.end();

    //kernels_->AssertClassification(event.queue(), outputs);
    event.emplace(outputs_token_, std::move(outputs));
    alpaka::wait(event.queue());
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "(ClassifierAot) E: " << event.id().event() << " OK - "
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
