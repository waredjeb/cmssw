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

    const float threshold = 0.6;
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
    auto const &gnn_input = event.get(gnn_input_token_);

    auto numEdges = gnn_output.view().metadata().size();
    auto edge_index_records = gnn_input.const_view<GNNEdgeIndexSoA>().records();
    GNNPostprocessingSoA::ConstView merged_view(
        gnn_output.const_view().records().score(), edge_index_records.in(), edge_index_records.out());

    TrackstersGNNPostprocessingSoAHostCollection gnn_post_host(numEdges, event.queue());
    gnn_post_host.deepCopy(merged_view, event.queue());
    alpaka::wait(event.queue());
    auto post_view = gnn_post_host.view();

    std::stringstream msg_stream;
    msg_stream << "MergedGNNTracksterProducer::produce [E: " << event.id().event() << "]";
   
    std::vector<ticl::Trackster> output(tracksters);
    std::vector<int> lookup(output.size());
    std::iota(lookup.begin(), lookup.end(), 0);
    std::array<int, 2> merge_idx;
    
    for (int i = 0; i < numEdges; i++) {
      if (post_view.score()[i] > threshold) {
        merge_idx[0] = post_view.out()[i];
        while (merge_idx[0] != lookup[merge_idx[0]]) {
          merge_idx[0] = lookup[merge_idx[0]];
        }

        merge_idx[1] = post_view.in()[i];
        while (merge_idx[1] != lookup[merge_idx[1]]) {
          merge_idx[1] = lookup[merge_idx[1]];
        }

        if (merge_idx[0] != merge_idx[1]) {
          output[merge_idx[0]].mergeTracksters(output[merge_idx[1]]);
          lookup[merge_idx[1]] = merge_idx[0];
        }
      }
    }

    for (int idx = lookup.size() - 1; idx >= 0; idx--) {
      if (lookup[idx] != idx) {
        output.erase(output.begin() + idx);
      }
    }

    output.shrink_to_fit();

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
