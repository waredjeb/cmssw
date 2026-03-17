#include <vector>
#include <stack>
#include <limits>
#include <Eigen/Core>

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/CLUE3DStateSoA.h"
#include "DataFormats/HGCalReco/interface/alpaka/CLUE3DStateDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersDeviceCollection.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "HeterogeneousCore/AlpakaCore/interface/MoveToDeviceCache.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/SynchronizingEDProducer.h"
#include "RecoHGCal/TICL/plugins/alpaka/CLUE3DKernel.h"
#include "RecoHGCal/TICL/plugins/alpaka/HGCalTiles.h"
#include "RecoHGCal/TICL/interface/alpaka/CLUE3DParamsSoA.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace ticl {

    namespace {
      using CLUE3DParamsCache =
          cms::alpakatools::MoveToDeviceCache<Device, PortableHostCollection<::ticl::CLUE3DParamsSoA>>;

      // Tile binning constants
      constexpr int kNEtaBins = 50;
      constexpr int kNPhiBins = 50;
      constexpr float kEtaMin = -3.5f;
      constexpr float kEtaMax = 3.5f;
      constexpr float kPhiMin = -M_PI;
      constexpr float kPhiMax = M_PI;

      inline int getEtaBin(float eta) {
        int bin = static_cast<int>((eta - kEtaMin) / (kEtaMax - kEtaMin) * kNEtaBins);
        return std::max(0, std::min(bin, kNEtaBins - 1));
      }

      inline int getPhiBin(float phi) {
        int bin = static_cast<int>((phi - kPhiMin) / (kPhiMax - kPhiMin) * kNPhiBins);
        return std::max(0, std::min(bin, kNPhiBins - 1));
      }
    }

    class CLUE3DProducer : public stream::SynchronizingEDProducer<edm::GlobalCache<CLUE3DParamsCache>> {
    public:
      static std::unique_ptr<CLUE3DParamsCache> initializeGlobalCache(edm::ParameterSet const& config) {
        PortableHostCollection<::ticl::CLUE3DParamsSoA> obj(cms::alpakatools::host(), ::ticl::kMaxAlgoIds);
        auto view = obj.view();

        // Load density calculation parameters
        auto criticalDensity = config.getParameter<std::vector<double>>("criticalDensity");
        auto criticalSelfDensity = config.getParameter<std::vector<double>>("criticalSelfDensity");
        auto densitySiblingLayers = config.getParameter<std::vector<int>>("densitySiblingLayers");
        auto densityEtaPhiDistanceSqr = config.getParameter<std::vector<double>>("densityEtaPhiDistanceSqr");
        auto densityXYDistanceSqr = config.getParameter<std::vector<double>>("densityXYDistanceSqr");
        auto kernelDensityFactor = config.getParameter<std::vector<double>>("kernelDensityFactor");

        // Load seed finding parameters
        auto criticalEtaPhiDistance = config.getParameter<std::vector<double>>("criticalEtaPhiDistance");
        auto criticalXYDistance = config.getParameter<std::vector<double>>("criticalXYDistance");
        auto criticalZDistanceLyr = config.getParameter<std::vector<int>>("criticalZDistanceLyr");
        auto outlierMultiplier = config.getParameter<std::vector<double>>("outlierMultiplier");
        auto minNumLayerCluster = config.getParameter<std::vector<int>>("minNumLayerCluster");

        // Populate arrays (assuming kMaxAlgoIds = 3)
        for (size_t i = 0; i < std::min(criticalDensity.size(), size_t(::ticl::kMaxAlgoIds)); ++i) {
          view.criticalDensity(i) = criticalDensity[i];
          view.criticalSelfDensity(i) = criticalSelfDensity[i];
          view.densitySiblingLayers(i) = densitySiblingLayers[i];
          view.densityEtaPhiDistanceSqr(i) = densityEtaPhiDistanceSqr[i];
          view.densityXYDistanceSqr(i) = densityXYDistanceSqr[i];
          view.kernelDensityFactor(i) = kernelDensityFactor[i];
          view.criticalEtaPhiDistance(i) = criticalEtaPhiDistance[i];
          view.criticalXYDistance(i) = criticalXYDistance[i];
          view.criticalZDistanceLyr(i) = criticalZDistanceLyr[i];
          view.outlierMultiplier(i) = outlierMultiplier[i];
          view.minNumLayerCluster(i) = minNumLayerCluster[i];
        }

        // Boolean flags
        view.densityOnSameLayer() = config.getParameter<bool>("densityOnSameLayer");
        view.nearestHigherOnSameLayer() = config.getParameter<bool>("nearestHigherOnSameLayer");
        view.useAbsoluteProjectiveScale() = config.getParameter<bool>("useAbsoluteProjectiveScale");
        view.useClusterDimensionXY() = config.getParameter<bool>("useClusterDimensionXY");
        view.rescaleDensityByZ() = config.getParameter<bool>("rescaleDensityByZ");

        // Geometry info (will be updated per-event, but set defaults here)
        view.lastLayerPerSide() = 50;  // typical HGCal, updated per-event
        view.nEtaBins() = config.getParameter<int>("nEtaBins");
        view.nPhiBins() = config.getParameter<int>("nPhiBins");
        view.etaMin() = -3.5f;
        view.etaMax() = 3.5f;
        view.phiMin() = -M_PI;
        view.phiMax() = M_PI;

        return std::make_unique<CLUE3DParamsCache>(std::move(obj));
      }

      CLUE3DProducer(edm::ParameterSet const& config, CLUE3DParamsCache const*)
          : SynchronizingEDProducer(config),
            inputLayerClusters_Token_{consumes(config.getParameter<edm::InputTag>("layerClusters"))},
            synchronise_(config.getParameter<bool>("synchronise")),
            algo_verbosity_(config.getParameter<int>("algo_verbosity")),
            minNumLayerCluster_(config.getParameter<std::vector<int>>("minNumLayerCluster")),
            nSeedsHost_(cms::alpakatools::make_host_buffer<int>()) {
        // Initialize RecHitTools (host-side)
        rhtools_ = std::make_unique<hgcal::RecHitTools>();

        // TODO: CaloGeometry access for RecHitTools
        // Currently commented out as CaloGeometry is host-only and not available
        // in device::EventSetup. Will need to be accessed differently.

        // TODO: Register output product properly
        // For now, this is a skeleton implementation without output
        // produces<std::vector<::ticl::Trackster>>();
      }

      void acquire(device::Event const& event, device::EventSetup const& setup) override {
        // NOTE: CaloGeometry is accessed from regular EventSetup (not device::EventSetup)
        // For now, we hardcode geometry parameters. In production, need proper setup.

        // Get input layer clusters
        auto const& layerClusters = event.get(inputLayerClusters_Token_);
        nClusters_ = layerClusters.size();

        if (nClusters_ == 0) {
          if (algo_verbosity_ > 0) {
            edm::LogInfo("CLUE3DProducer") << "No layer clusters, skipping";
          }
          return;
        }

        // For geometry, we need lastLayerPerSide
        // TODO: Get from rhtools_ which needs CaloGeometry
        lastLayerPerSide_ = 50;  // Hardcoded for now
        int nLayers = 2 * lastLayerPerSide_;

        // 1. ALLOCATE DEVICE STATE
        clue3dState_.emplace(event.queue(), nClusters_);

        // 2. PREPARE HOST DATA
        // Create host mirror of state for initialization
        PortableHostCollection<CLUE3DStateSoA> hostState(cms::alpakatools::host(), nClusters_);
        auto hostView = hostState.view();

        layerIndices_.resize(nClusters_);

        // 3. POPULATE STATE FROM LAYER CLUSTERS
        for (int i = 0; i < nClusters_; ++i) {
          auto const& lc = layerClusters[i];

          // Determine layer (simplified - normally use rhtools_)
          // TODO: Get from detector ID via rhtools_
          int layer = 0;  // Placeholder
          layerIndices_[i] = layer;

          // Calculate radius (simplified version without rhtools for now)
          float radius = 1.0f;  // TODO: Use calculateClusterRadius(lc) when geometry is available

          // Fill host state
          hostView.x(i) = lc.x();
          hostView.y(i) = lc.y();
          hostView.z(i) = lc.z();
          hostView.eta(i) = lc.eta();
          hostView.phi(i) = lc.phi();
          hostView.r_over_absz(i) = std::sqrt(lc.x() * lc.x() + lc.y() * lc.y()) / std::abs(lc.z());
          hostView.radius(i) = radius;
          hostView.energy(i) = lc.energy();
          hostView.cells(i) = lc.hitsAndFractions().size();
          hostView.layer(i) = layer;
          hostView.algoId(i) = lc.algo() - reco::CaloCluster::hgcal_em;  // 0, 1, or 2
          hostView.isSilicon(i) = 1;  // TODO: Get from detector ID
          hostView.layerClusterOriginalIdx(i) = i;

          // Initialize clustering state
          hostView.rho(i) = 0.0f;
          hostView.z_extension(i) = 0.0f;
          hostView.delta_dist(i) = std::numeric_limits<float>::max();
          hostView.delta_layer(i) = std::numeric_limits<int>::max();
          hostView.nearestHigher_layer(i) = -1;
          hostView.nearestHigher_idx(i) = -1;
          hostView.clusterIndex(i) = -1;
          hostView.isSeed(i) = 0;
          hostView.isOutlier(i) = 0;
        }

        // 4. TRANSFER STATE TO DEVICE
        alpaka::memcpy(event.queue(), clue3dState_->buffer(), hostState.buffer());

        // 5. CREATE AND FILL TILES ON DEVICE
        tiles_.emplace(event.queue(), nClusters_, nLayers, kNEtaBins, kNPhiBins);
        tiles_->fill(event.queue(), *clue3dState_, nClusters_);

        // 6. PREPARE LAYER Z POSITIONS (TODO)
        // For now, we'll skip layer Z positions and handle in kernels differently
        float* layersZ_ptr = nullptr;  // Placeholder

        // 7. PREPARE SEED COUNTER
        auto nSeeds_d = cms::alpakatools::make_device_buffer<int>(event.queue());
        alpaka::memset(event.queue(), nSeeds_d, 0);

        // 8. GET PARAMETERS
        auto const& params = globalCache()->get(event.queue());

        // 9. LAUNCH KERNELS
        CLUE3DKernel kernel;

        const auto& tilesView = tiles_->view();

        kernel.calculateLocalDensity(
            event.queue(), params.const_view(), tilesView, layersZ_ptr, *clue3dState_, nClusters_);

        kernel.calculateDistanceToHigher(
            event.queue(), params.const_view(), tilesView, *clue3dState_, nClusters_);

        kernel.findAndAssignSeeds(event.queue(), params.const_view(), *clue3dState_, nClusters_, alpaka::getPtrNative(nSeeds_d));

        // 10. TRANSFER SEED COUNT BACK
        alpaka::memcpy(event.queue(), nSeedsHost_, nSeeds_d);
        alpaka::wait(event.queue());  // Wait for seed count

        nSeeds_ = *alpaka::getPtrNative(nSeedsHost_);

        if (algo_verbosity_ > 0) {
          edm::LogInfo("CLUE3DProducer") << "Found " << nSeeds_ << " seeds from " << nClusters_ << " clusters";
        }
      }

      void produce(device::Event& event, device::EventSetup const& setup) override {
        // Handle empty events
        if (!clue3dState_ || nClusters_ == 0) {
          // Produce empty collection
          // TODO: Register and emit output when output token is added
          if (synchronise_)
            alpaka::wait(event.queue());
          return;
        }

        // 1. TRANSFER STATE BACK TO HOST
        PortableHostCollection<CLUE3DStateSoA> hostState(cms::alpakatools::host(), nClusters_);
        alpaka::memcpy(event.queue(), hostState.buffer(), clue3dState_->buffer());
        alpaka::wait(event.queue());

        auto stateView = hostState.view();

        // 2. BUILD FOLLOWER GRAPH ON HOST
        std::vector<std::vector<int>> followers(nClusters_);
        std::vector<int> tracksterSeedAlgoId;

        for (int i = 0; i < nClusters_; ++i) {
          if (!stateView.isSeed(i) && !stateView.isOutlier(i)) {
            int higher = stateView.nearestHigher_idx(i);
            if (higher >= 0 && higher < nClusters_) {
              followers[higher].push_back(i);
            }
          }
          if (stateView.isSeed(i)) {
            tracksterSeedAlgoId.push_back(stateView.algoId(i));
          }
        }

        // 3. PROPAGATE CLUSTER INDICES VIA DFS
        std::stack<int> stack;
        for (int i = 0; i < nClusters_; ++i) {
          if (stateView.isSeed(i)) {
            stack.push(i);
          }
        }

        while (!stack.empty()) {
          int idx = stack.top();
          stack.pop();
          int clusterIdx = stateView.clusterIndex(idx);

          for (int follower : followers[idx]) {
            stateView.clusterIndex(follower) = clusterIdx;
            stack.push(follower);
          }
        }

        // 4. BUILD TRACKSTER VERTICES AND EDGES
        std::vector<std::vector<unsigned int>> tracksterVertices(nSeeds_);
        std::vector<std::vector<std::array<unsigned int, 2>>> tracksterEdges(nSeeds_);

        for (int i = 0; i < nClusters_; ++i) {
          int tIdx = stateView.clusterIndex(i);
          if (tIdx >= 0 && tIdx < nSeeds_) {
            int origIdx = stateView.layerClusterOriginalIdx(i);
            tracksterVertices[tIdx].push_back(origIdx);

            // Add edges to followers
            for (int followerLocal : followers[i]) {
              int followerOrig = stateView.layerClusterOriginalIdx(followerLocal);
              tracksterEdges[tIdx].push_back(
                  {{static_cast<unsigned int>(origIdx), static_cast<unsigned int>(followerOrig)}});
            }
          }
        }

        // 5. CREATE TRACKSTER OBJECTS
        std::vector<::ticl::Trackster> tracksters;
        for (int t = 0; t < nSeeds_; ++t) {
          int algoId = tracksterSeedAlgoId[t];
          if (tracksterVertices[t].size() >= static_cast<size_t>(minNumLayerCluster_[algoId])) {
            ::ticl::Trackster trackster;
            trackster.vertices() = tracksterVertices[t];
            trackster.vertex_multiplicity().resize(tracksterVertices[t].size(), 1);
            trackster.edges() = tracksterEdges[t];

            tracksters.push_back(trackster);
          }
        }

        if (algo_verbosity_ > 0) {
          edm::LogInfo("CLUE3DProducer") << "Created " << tracksters.size() << " tracksters";
        }

        // 6. COMPUTE PCA (TODO)
        // This requires access to layer clusters
        // For now, tracksters have vertices and edges but no PCA properties

        // 7. OUTPUT (TODO)
        // event.emplace(outputToken, std::move(tracksters));

        if (synchronise_)
          alpaka::wait(event.queue());
      }

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
        edm::ParameterSetDescription desc;

        desc.add<edm::InputTag>("layerClusters", edm::InputTag(""));

        desc.add<std::vector<double>>("criticalDensity", {4., 4., 4.})->setComment("in GeV");
        desc.add<std::vector<double>>("criticalSelfDensity", {0.15, 0.15, 0.15})
            ->setComment("Minimum ratio of self_energy/local_density to become a seed");
        desc.add<std::vector<int>>("densitySiblingLayers", {3, 3, 3})
            ->setComment("inclusive, layers to consider while computing local density and searching for nearestHigher");
        desc.add<std::vector<double>>("densityEtaPhiDistanceSqr", {0.0008, 0.0008, 0.0008})
            ->setComment("in eta,phi space, distance to consider for local density");
        desc.add<std::vector<double>>("densityXYDistanceSqr", {3.24, 3.24, 3.24})
            ->setComment("in cm, distance on the transverse plane to consider for local density");
        desc.add<std::vector<double>>("kernelDensityFactor", {0.2, 0.2, 0.2})
            ->setComment("Kernel factor to be applied to other LC while computing the local density");
        desc.add<bool>("densityOnSameLayer", false);
        desc.add<bool>("nearestHigherOnSameLayer", false)
            ->setComment("Allow the nearestHigher to be located on the same layer");
        desc.add<bool>("useAbsoluteProjectiveScale", true)
            ->setComment("Express all cuts in terms of r/z*z_0{,phi} projective variables");
        desc.add<bool>("useClusterDimensionXY", false)
            ->setComment(
                "If true use the estimated cluster radius to determine compatibility while computing local density");
        desc.add<bool>("rescaleDensityByZ", false)
            ->setComment("Rescale local density by the Z volume explored");

        desc.add<std::vector<double>>("criticalEtaPhiDistance", {0.025, 0.025, 0.025})
            ->setComment("Minimal distance in eta,phi from nearestHigher to become a seed");
        desc.add<std::vector<double>>("criticalXYDistance", {1.8, 1.8, 1.8})
            ->setComment("Minimal distance in cm on XY plane from nearestHigher to become a seed");
        desc.add<std::vector<int>>("criticalZDistanceLyr", {5, 5, 5})
            ->setComment("Minimal distance in layers along Z from nearestHigher to become a seed");
        desc.add<std::vector<double>>("outlierMultiplier", {2., 2., 2.})
            ->setComment("Multiplier for transverse distance to classify as outlier");
        desc.add<std::vector<int>>("minNumLayerCluster", {2, 2, 2})->setComment("Minimum layer clusters per trackster");

        desc.add<int>("nEtaBins", 50)->setComment("Number of eta bins for tiles");
        desc.add<int>("nPhiBins", 50)->setComment("Number of phi bins for tiles");

        desc.add<int>("algo_verbosity", 0);
        desc.add<bool>("synchronise", false);

        descriptions.addWithDefaultLabel(desc);
      }

      static void globalEndJob(CLUE3DParamsCache*) {}

    private:
      const edm::EDGetTokenT<std::vector<reco::CaloCluster>> inputLayerClusters_Token_;
      const bool synchronise_;
      const int algo_verbosity_;
      const std::vector<int> minNumLayerCluster_;

      std::unique_ptr<hgcal::RecHitTools> rhtools_;

      // Device collections (created per-event in acquire)
      std::optional<CLUE3DStateDeviceCollection> clue3dState_;
      std::optional<HGCalTiles> tiles_;

      // Per-event state
      int nSeeds_ = 0;
      int nClusters_ = 0;
      int lastLayerPerSide_ = 0;
      std::vector<int> layerIndices_;

      // Host buffers
      cms::alpakatools::host_buffer<int> nSeedsHost_;

      // Helper method
      float calculateClusterRadius(const reco::CaloCluster& lc) const {
        float sum_x = 0.f, sum_y = 0.f, sum_sqr_x = 0.f, sum_sqr_y = 0.f;
        float ref_x = lc.x(), ref_y = lc.y();
        float invClsize = 1.f / lc.hitsAndFractions().size();

        for (auto const& hf : lc.hitsAndFractions()) {
          auto const& point = rhtools_->getPosition(hf.first);
          sum_x += point.x() - ref_x;
          sum_sqr_x += (point.x() - ref_x) * (point.x() - ref_x);
          sum_y += point.y() - ref_y;
          sum_sqr_y += (point.y() - ref_y) * (point.y() - ref_y);
        }

        float radius_x = std::sqrt((sum_sqr_x - (sum_x * sum_x) * invClsize) * invClsize);
        float radius_y = std::sqrt((sum_sqr_y - (sum_y * sum_y) * invClsize) * invClsize);

        // Handle single-cell clusters
        if (invClsize == 1.f) {
          auto detId = lc.hitsAndFractions()[0].first;
          if (rhtools_->isSilicon(detId)) {
            radius_x = radius_y = rhtools_->getRadiusToSide(detId);
          } else {
            auto const& point = rhtools_->getPosition(detId);
            auto const& eta_phi_window = rhtools_->getScintDEtaDPhi(detId);
            radius_x = radius_y = point.perp() * eta_phi_window.second;
          }
        }

        return radius_x + radius_y;
      }

      // TODO: Add CaloGeometry token (needs host-side access)
      // TODO: Add output token when implementing full data flow
    };

  }  // namespace ticl
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(ticl::CLUE3DProducer);
