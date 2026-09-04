#include <algorithm>
#include <cstdint>
#include <optional>
#include <vector>

#include "DataFormats/CaloRecHit/interface/alpaka/CaloClusterDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/SynchronizingEDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

#include "TracksterCLUEsteringAlgoWrapper.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // Run CLUEstering on CaloClusterDeviceCollection, alread on device from Layer cluster step.
  // Returns the tracster index of each layer cluster in the merged host collection
  // Trackster properties are filled on host by ticl::PatternRecognitionbyCLUEstering
  class TrackstersCLUEsteringProducer : public stream::SynchronizingEDProducer<> {
  public:
    TrackstersCLUEsteringProducer(edm::ParameterSet const& config)
        : SynchronizingEDProducer(config),
          maskToken_{consumes<std::vector<float>>(config.getParameter<edm::InputTag>("filtered_mask"))},
          assignmentToken_{produces()},
          parameters_{static_cast<float>(config.getParameter<double>("dc")),
                      static_cast<float>(config.getParameter<double>("rhoc")),
                      static_cast<float>(config.getParameter<double>("outlierDistance")),
                      static_cast<float>(config.getParameter<double>("seedingDistance")),
                      [&config] {
                        auto const& v = config.getParameter<std::vector<double>>("sigmaT");
                        return std::vector<float>(v.begin(), v.end());
                      }(),
                      static_cast<float>(config.getParameter<double>("layerScale")),
                      static_cast<float>(config.getParameter<double>("rhocEtaExponent")),
                      static_cast<float>(config.getParameter<double>("rhocPivotRadius"))} {
      // CLUEstering requires seedingDistance <= outlierDistance
      if (parameters_.seedingDistance > parameters_.outlierDistance) {
        throw edm::Exception(edm::errors::Configuration)
            << "TrackstersCLUEsteringProducer: seedingDistance (" << parameters_.seedingDistance
            << ") must not exceed outlierDistance (" << parameters_.outlierDistance << ").";
      }
      // The order of these tags must match the one used by MergeClusterProducer to build the
      // merged host collection, because the indices produced here index into it.
      for (const auto& tag : config.getParameter<std::vector<edm::InputTag>>("layerClusters")) {
        layerClustersTokens_.emplace_back(consumes(tag));
      }
      if (parameters_.sigmaT.size() != layerClustersTokens_.size()) {
        throw edm::Exception(edm::errors::Configuration)
            << "TrackstersCLUEsteringProducer: 'sigmaT' has " << parameters_.sigmaT.size()
            << " entries but 'layerClusters' lists " << layerClustersTokens_.size()
            << " collections; one sigmaT per collection is required, in the same order.";
      }
      for (auto s : parameters_.sigmaT) {
        if (s <= 0.f) {
          throw edm::Exception(edm::errors::Configuration)
              << "TrackstersCLUEsteringProducer: every sigmaT entry must be positive.";
        }
      }
    }

    ~TrackstersCLUEsteringProducer() override = default;

    void acquire(device::Event const& iEvent, device::EventSetup const& iSetup) override {
      std::vector<::reco::CaloClusterSoAConstView> views;
      std::vector<uint32_t> offsets;
      views.reserve(layerClustersTokens_.size());
      offsets.reserve(layerClustersTokens_.size());

      uint32_t nTotal = 0;
      for (auto const& token : layerClustersTokens_) {
        auto const& collection = iEvent.get(token);
        views.emplace_back(collection.view());
        offsets.emplace_back(nTotal);
        // CaloClusterSoA is a multi-block SoA, so metadata().size() is one size per block.
        // All blocks are allocated with the number of clusters, so block 0 is that number.
        nTotal += static_cast<uint32_t>(collection.view().metadata().size()[0]);
      }

      auto const& mask = iEvent.get(maskToken_);
      if (mask.size() != nTotal) {
        throw edm::Exception(edm::errors::Configuration)
            << "TrackstersCLUEsteringProducer: the filtered layer-cluster mask has " << mask.size()
            << " entries but the layer-cluster device collections hold " << nTotal
            << " clusters in total. The 'layerClusters' tags must list exactly the collections "
               "that MergeClusterProducer merges, in the same order.";
      }

      nTotal_ = nTotal;
      // Nothing to cluster: the buffers below stay unallocated and produce() publishes an empty
      // assignment.
      if (nTotal == 0) {
        return;
      }

      // CLUEstering needs point count to size its tiles, calculating here avoids a device to
      // host transfer. There may be a better solution
      const auto nSurviving =
          static_cast<uint32_t>(std::count_if(mask.begin(), mask.end(), [](float m) { return m > 0.f; }));

      auto& queue = iEvent.queue();
      deviceMask_ = cms::alpakatools::make_device_buffer<float[]>(queue, nTotal);
      deviceAssignment_ = cms::alpakatools::make_device_buffer<int32_t[]>(queue, nTotal);
      hostAssignment_ = cms::alpakatools::make_host_buffer<int32_t[]>(queue, nTotal);

      auto hostMask = cms::alpakatools::make_host_view<const float>(mask.data(), nTotal);
      alpaka::memcpy(queue, *deviceMask_, hostMask);

      algo_.run(queue,
                std::span<const ::reco::CaloClusterSoAConstView>(views),
                std::span<const uint32_t>(offsets),
                deviceMask_->data(),
                nTotal,
                nSurviving,
                parameters_,
                deviceAssignment_->data());

      alpaka::memcpy(queue, *hostAssignment_, *deviceAssignment_);
    }

    void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {
      // acquire() returns early on an empty event, so the buffers are only engaged when nTotal_ > 0.
      std::vector<int32_t> assignment;
      if (nTotal_ > 0) {
        assignment.assign(hostAssignment_->data(), hostAssignment_->data() + nTotal_);
      }
      iEvent.emplace(assignmentToken_, std::move(assignment));
      deviceMask_.reset();
      deviceAssignment_.reset();
      hostAssignment_.reset();
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<std::vector<edm::InputTag>>("layerClusters",
                                           {edm::InputTag("hgcalSoALayerClustersEE"),
                                            edm::InputTag("hgcalSoALayerClustersHSi"),
                                            edm::InputTag("hgcalSoALayerClustersHSci")})
          ->setComment(
              "Per-detector layer-cluster SoA on device, in the same order as the 'layerClusters' "
              "tags of MergeClusterProducer.");
      desc.add<edm::InputTag>("filtered_mask", edm::InputTag("filteredLayerClustersCLUE3DHigh", "CLUE3DHigh"));
      desc.add<double>("dc", 1.)->setComment(
          "Density radius. The clustering runs in (x/|z|/sigmaT, y/|z|/sigmaT, layer/layerScale) with the "
          "cylinder metric max(|dT|, |dl|), so dc, seedingDistance and outlierDistance are in units of "
          "sigmaT transversely and layerScale longitudinally.");
      desc.add<double>("rhoc", 0.8)->setComment("Minimum energy density for a point to seed a trackster.");
      desc.add<double>("outlierDistance", 3.6)->setComment("Minimum distance between clusters.");
      desc.add<double>("seedingDistance", 2.8)->setComment("Distance threshold for seed points.");
      desc.add<std::vector<double>>("sigmaT", {0.003, 0.006, 0.012})
          ->setComment(
              "Angular transverse scale [rad], one entry per 'layerClusters' collection (EE, HSi, HSci "
              "with the default tags), in the same order.");
      desc.add<double>("layerScale", 4.)->setComment("Longitudinal scale in layers.");
      desc.add<double>("rhocEtaExponent", 2.)
          ->setComment(
              "Exponent alpha of the eta-dependent seed threshold rhoc(point) = rhoc * (rhocPivotRadius/r)^alpha, "
              "with r the point's gnomonic radius sqrt(x^2+y^2)/|z|. 0 keeps rhoc flat.");
      desc.add<double>("rhocPivotRadius", 0.42)
          ->setComment(
              "Gnomonic radius where the eta-dependent seed threshold equals rhoc (0.42 corresponds to eta ~ 1.6).");
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    std::vector<device::EDGetToken<reco::CaloClusterDeviceCollection>> layerClustersTokens_;
    edm::EDGetTokenT<std::vector<float>> const maskToken_;
    edm::EDPutTokenT<std::vector<int32_t>> const assignmentToken_;

    TracksterCLUEsteringAlgoWrapper algo_;
    TracksterCLUEsteringAlgoWrapper::Parameters const parameters_;

    std::optional<cms::alpakatools::device_buffer<Device, float[]>> deviceMask_;
    std::optional<cms::alpakatools::device_buffer<Device, int32_t[]>> deviceAssignment_;
    std::optional<cms::alpakatools::host_buffer<int32_t[]>> hostAssignment_;
    uint32_t nTotal_ = 0;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(TrackstersCLUEsteringProducer);
