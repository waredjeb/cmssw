#include "DataFormats/HGCalReco/interface/alpaka/TracksterSoADeviceCollection.h"
#include "DataFormats/HGCalReco/interface/TracksterSoAHostCollection.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLGraph.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "PhysicsTools/PyTorch/interface/Nvtx.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class MergedGNNTracksterProducer : public stream::EDProducer<> {
  public:
    MergedGNNTracksterProducer(const edm::ParameterSet &params);

    void produce(device::Event &event, const device::EventSetup &event_setup) override;
    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  private:
    const edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_token_;
    const device::EDGetToken<TrackstersSoADeviceCollection> gnn_input_token_;
    const device::EDGetToken<TrackstersGNNOutputSoADeviceCollection> gnn_output_token_;
    const edm::EDPutTokenT<std::vector<ticl::Trackster>> merged_tracksters_token_; /**< Token to store output data. */

    const int threshold = 0.6;
  };

  MergedGNNTracksterProducer::MergedGNNTracksterProducer(edm::ParameterSet const &params)
      : EDProducer<>(params),
        tracksters_token_(consumes<std::vector<ticl::Trackster>>(params.getParameter<edm::InputTag>("tracksters"))),
        gnn_input_token_{consumes(params.getParameter<edm::InputTag>("gnnInput"))},
        gnn_output_token_(consumes(params.getParameter<edm::InputTag>("gnnOutput"))),
        merged_tracksters_token_{produces()} {}

  void MergedGNNTracksterProducer::produce(device::Event &event, const device::EventSetup &event_setup) {
    auto const &tracksters = event.get(tracksters_token_);
    auto const &gnn_output = event.get(gnn_output_token_);
    auto const &gnn_input = const_cast<TrackstersSoADeviceCollection &>(event.get(gnn_input_token_));

    const int numNodes = gnn_input.const_view<GNNNodeSoA>().metadata().size();
    const int numEdges = gnn_input.const_view<GNNEdgeSoA>().metadata().size();
    std::array<int, 3> const sizes{{numNodes, numEdges, numEdges}};

    auto gnn_input_host = TrackstersSoAHostCollection(sizes, event.queue());
    auto gnn_output_host = TrackstersGNNOutputSoAHostCollection(numEdges, event.queue());

    alpaka::memcpy(event.queue(), gnn_output_host.buffer(), gnn_output.buffer());
    alpaka::memcpy(event.queue(), gnn_input_host.buffer(), gnn_input.buffer());
    alpaka::wait(event.queue());

    std::stringstream msg_stream;
    msg_stream << "MergedGNNTracksterProducer::produce [E: " << event.id().event() << "]";

    std::vector<ticl::Node> nodes;
    for (size_t i = 0; i < tracksters.size(); i++) {
      nodes.emplace_back(i);
    }

    TICLGraph ticlGraph(nodes);

    auto outputs_view = gnn_output_host.view();
    auto inputs_view = gnn_input_host.view<GNNEdgeIndexSoA>();

    for (int i = 0; i < outputs_view.metadata().size(); i++) {
      if (outputs_view.score()[i] > threshold) {
        ticlGraph.adaptNode(inputs_view.out()[i]).addOuterNeighbour(inputs_view.in()[i]);
        ticlGraph.adaptNode(inputs_view.in()[i]).addInnerNeighbour(inputs_view.out()[i]);
      }
    }

    ticlGraph.findRootNodes();

    std::vector<ticl::Trackster> output;
    std::vector<ticl::Trackster> tmp;
    ticl::Trackster trackster;
    auto components = ticlGraph.findSubComponents();
    for (auto comp : components) {
      if (comp.size() < 1) {
        continue;
      }
      trackster = tracksters[comp.back()];

      if (comp.size() > 1) {
        comp.pop_back();
        for (auto node : comp) {
          tmp.push_back(tracksters[node]);
        }
        trackster.mergeTracksters(tmp);
        tmp.clear();
      }
      output.push_back(trackster);
    }

    std::cout << "(MergedGNNTracksterProducer) Number of Trackster: " << output.size() << std::endl;
    event.emplace(merged_tracksters_token_, std::move(output));
  }

  /**
   * @brief Describes the allowed configuration parameters for this module.
   * @param descriptions Configuration description object to populate.
   */
  void MergedGNNTracksterProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracksters", edm::InputTag("ticlTrackstersCLUE3DHigh"));
    desc.add<edm::InputTag>("gnnInput");
    desc.add<edm::InputTag>("gnnOutput");
    descriptions.addWithDefaultLabel(desc);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(MergedGNNTracksterProducer);
