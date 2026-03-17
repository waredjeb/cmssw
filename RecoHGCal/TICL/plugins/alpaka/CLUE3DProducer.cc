#include <vector>
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
#include "DataFormats/HGCalReco/interface/alpaka/CLUE3DStateDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalTilesDeviceCollection.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "HeterogeneousCore/AlpakaCore/interface/MoveToDeviceCache.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/SynchronizingEDProducer.h"
#include "RecoHGCal/TICL/plugins/alpaka/CLUE3DKernel.h"
#include "RecoHGCal/TICL/interface/alpaka/CLUE3DParamsSoA.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace ticl {

    namespace {
      using CLUE3DParamsCache =
          cms::alpakatools::MoveToDeviceCache<Device, PortableHostCollection<::ticl::CLUE3DParamsSoA>>;
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
            algo_verbosity_(config.getParameter<int>("algo_verbosity")) {
        // Initialize RecHitTools (host-side)
        rhtools_ = std::make_unique<hgcal::RecHitTools>();

        // TODO: CaloGeometry access for RecHitTools
        // Currently commented out as CaloGeometry is host-only and not available
        // in device::EventSetup. Will need to be accessed differently.

        // TODO: Register output product properly
        // For now, this is a skeleton implementation without output
      }

      void acquire(device::Event const& event, device::EventSetup const& setup) override {
        // TODO: For first version, this acquire method will prepare data on host
        // In a full implementation, we'd accept device collections directly if available
        // For now: stub that shows the structure

        // TODO: Get geometry (needs host-side access, not available in device::EventSetup)
        // const CaloGeometry& geom = ...;  // Need to get from regular EventSetup
        // rhtools_->setGeometry(geom);
        // int lastLayerPerSide = rhtools_->lastLayer(false);

        // TODO: Get input layer clusters (assuming host collection for now)
        // In full implementation, this would be a device collection or we'd transfer
        // auto const& layerClusters = event.get(inputLayerClusters_Token_);

        // TODO: Build tiles structure on host, transfer to device
        // TODO: Prepare CLUE3D state collection on device
        // TODO: Launch kernels via CLUE3DKernel

        // For now, just log that acquire was called
        if (algo_verbosity_ > 0) {
          edm::LogInfo("CLUE3DProducer") << "acquire() called (skeleton)";
        }
      }

      void produce(device::Event& event, device::EventSetup const& setup) override {
        // TODO: Implement full produce logic:
        // 1. Transfer results from device to host
        // 2. Run host-side graph propagation (cluster index assignment)
        // 3. Build Trackster objects
        // 4. Apply filters (min layer cluster count, PID cuts)
        // 5. Register and produce output to event
        //
        // NOTE: Output product registration needs to be added once the full
        // data flow is implemented. The Alpaka producer interface requires
        // careful handling of device vs host collections.

        if (algo_verbosity_ > 0) {
          edm::LogInfo("CLUE3DProducer") << "produce() called (skeleton - no output yet)";
        }

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

      std::unique_ptr<hgcal::RecHitTools> rhtools_;

      // Device collections (created per-event in acquire)
      std::optional<CLUE3DStateDeviceCollection> clue3dState_;
      std::optional<HGCalTilesDeviceCollection> tiles_;

      // TODO: Add CaloGeometry token (needs host-side access)
      // TODO: Add output token when implementing full data flow
    };

  }  // namespace ticl
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(ticl::CLUE3DProducer);
