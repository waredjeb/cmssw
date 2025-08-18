#include <memory>
#include "DataFormats/HGCalReco/interface/alpaka/TracksterSoADeviceCollection.h"
#include "DataFormats/HGCalReco/interface/TracksterSoAHostCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/HGCalReco/interface/TICLGraph.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/HGCalReco/interface/Common.h"
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
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "TrackstersPCAAlpaka.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class MergedGNNTracksterProducer : public stream::EDProducer<> {
  public:
    MergedGNNTracksterProducer(const edm::ParameterSet &params);

    void produce(device::Event &event, const device::EventSetup &event_setup) override;
    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
    void beginRun(edm::Run const &iEvent, edm::EventSetup const &es) override;

  private:
    const edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_token_;
    const device::EDGetToken<TrackstersSoADeviceCollection> gnn_input_token_;
    const device::EDGetToken<TrackstersGNNOutputSoADeviceCollection> gnn_output_token_;
    const edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;
    const edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> clustersTime_token_;
    const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometry_token_;
    const std::string detector_;
    const std::string propName_;

    const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfield_token_;
    const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagator_token_;
    const edm::EDPutTokenT<std::vector<ticl::Trackster>> merged_tracksters_token_; /**< Token to store output data. */
    const edm::EDPutTokenT<std::vector<std::vector<unsigned int>>> linked_merged_trackstersId_token_; /**< Token to store output data. */
    const HGCalDDDConstants *hgcons_;
    hgcal::RecHitTools rhtools_;
    edm::ESGetToken<HGCalDDDConstants, IdealGeometryRecord> hdc_token_;
  };

  MergedGNNTracksterProducer::MergedGNNTracksterProducer(edm::ParameterSet const &params)
      : EDProducer<>(params),
        tracksters_token_(consumes<std::vector<ticl::Trackster>>(params.getParameter<edm::InputTag>("tracksters"))),
        gnn_input_token_{consumes(params.getParameter<edm::InputTag>("gnnInput"))},
        gnn_output_token_(consumes(params.getParameter<edm::InputTag>("gnnOutput"))),
        clusters_token_(consumes<std::vector<reco::CaloCluster>>(params.getParameter<edm::InputTag>("layer_clusters"))),
        clustersTime_token_(
            consumes<edm::ValueMap<std::pair<float, float>>>(params.getParameter<edm::InputTag>("layer_clustersTime"))),
        geometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord, edm::Transition::BeginRun>()),
        detector_(params.getParameter<std::string>("detector")),
        propName_(params.getParameter<std::string>("propagator")),
        bfield_token_(esConsumes<MagneticField, IdealMagneticFieldRecord, edm::Transition::BeginRun>()),
        propagator_token_(esConsumes<Propagator, TrackingComponentsRecord, edm::Transition::BeginRun>(
            edm::ESInputTag("", propName_))),
        merged_tracksters_token_{produces()},
        linked_merged_trackstersId_token_{produces()}
  {
    std::string detectorName_ = (detector_ == "HFNose") ? "HGCalHFNoseSensitive" : "HGCalEESensitive";
    hdc_token_ = esConsumes<HGCalDDDConstants, IdealGeometryRecord, edm::Transition::BeginRun>(
        edm::ESInputTag("", detectorName_));
  }
  void MergedGNNTracksterProducer::beginRun(edm::Run const &iEvent, edm::EventSetup const &es) {
    edm::ESHandle<HGCalDDDConstants> hdc = es.getHandle(hdc_token_);
    hgcons_ = hdc.product();

    edm::ESHandle<CaloGeometry> geom = es.getHandle(geometry_token_);
    rhtools_.setGeometry(*geom);

    edm::ESHandle<MagneticField> bfield = es.getHandle(bfield_token_);
    edm::ESHandle<Propagator> propagator = es.getHandle(propagator_token_);
  };

  void MergedGNNTracksterProducer::produce(device::Event &event, const device::EventSetup &event_setup) {
    auto const &tracksters = event.get(tracksters_token_);
    auto const &gnn_output = event.get(gnn_output_token_);
    auto const &gnn_input = event.get(gnn_input_token_);
    const auto &layerClusters = event.get(clusters_token_);
    const auto &layerClustersTimes = event.get(clustersTime_token_);

    auto numEdges = gnn_output.view().metadata().size();
    auto edge_index_records = gnn_input.const_view<GNNEdgeIndexSoA>().records();
    GNNPostprocessingSoA::ConstView merged_view(
        gnn_output.const_view().records().score(), edge_index_records.in(), edge_index_records.out());

    TrackstersGNNPostprocessingSoAHostCollection gnn_post_host(numEdges, event.queue());
    gnn_post_host.deepCopy(merged_view, event.queue());
    alpaka::wait(event.queue());
    auto post_view = gnn_post_host.view();

    std::stringstream msg_stream;
    std::vector<ticl::Trackster> output(tracksters);
    std::vector<std::vector<unsigned int>> linkedTrackstersOutput(tracksters.size());
    std::vector<int> lookup(output.size());
    std::iota(lookup.begin(), lookup.end(), 0);
    std::array<int, 2> merge_idx;
    for (int i = 0; i < numEdges; ++i) {
      if (post_view.score()[i] > 0.99) {
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
    int nextIdx = 0;
    for (int i = 0; i < static_cast<int>(lookup.size()); ++i) {
      if (lookup[i] == i) {
        linkedTrackstersOutput[i].push_back(static_cast<unsigned int>(nextIdx));
        ++nextIdx;
      }
    }
    for (int idx = static_cast<int>(lookup.size()) - 1; idx >= 0; --idx) {
      if (lookup[idx] != idx) {
        output.erase(output.begin() + idx);
        linkedTrackstersOutput.erase(linkedTrackstersOutput.begin() + idx);
      }
    }

    ticlAlpaka::assignPCAtoTracksters(output,
                                layerClusters,
                                layerClustersTimes,
                                rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z(),
                                rhtools_,
                                true);
    event.emplace(merged_tracksters_token_, std::move(output));
    event.emplace(linked_merged_trackstersId_token_, std::move(linkedTrackstersOutput));
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
    desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalMergeLayerClusters"));
    desc.add<edm::InputTag>("layer_clustersTime", edm::InputTag("hgcalMergeLayerClusters", "timeLayerCluster"));
    desc.add<std::string>("detector", "HGCAL");
    desc.add<std::string>("propagator", "PropagatorWithMaterial");
    descriptions.addWithDefaultLabel(desc);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(MergedGNNTracksterProducer);
