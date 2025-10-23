#include <alpaka/alpaka.hpp>
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "DataFormats/HGCalReco/interface/HGCalSoAClusters.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "CLUEstering/CLUEstering.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class HeterogeneousTracksterProducer : public stream::EDProducer<> {
  public:
    HeterogeneousTracksterProducer(edm::ParameterSet const& config)
        : EDProducer(config),
          deviceTokenSoAClusters_{consumes(config.getParameter<edm::InputTag>("layerClusters"))},
          legacyTrackstersToken_{produces()},
          rho_(config.getParameter<double>("rho_c")) {
      auto dc_vec = config.getParameter<std::vector<double>>("dc");
      auto dm_vec = config.getParameter<std::vector<double>>("dm");

      if (dc_vec.size() != 3 || dm_vec.size() != 3) {
        throw cms::Exception("Configuration") << "Parameters 'dc' and 'dm' must each have exactly 3 elements.";
      }

      for (size_t i = 0; i < 3; ++i) {
        dc_[i] = static_cast<float>(dc_vec[i]);
        dm_[i] = static_cast<float>(dm_vec[i]);
      }
    }
    ~HeterogeneousTracksterProducer() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalRecHitsLayerClustersSoA"));
      desc.add<double>("rho_c", 0.6);
      desc.add<std::vector<double>>("dc", {2., 2., 2});
      desc.add<std::vector<double>>("dm", {1.8, 1.8, 2});
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {
      const auto& lc = iEvent.get(deviceTokenSoAClusters_);
      auto& queue = iEvent.queue();
      auto x = const_cast<float*>(lc.view().x().data());
      auto y = const_cast<float*>(lc.view().y().data());
      auto z = const_cast<float*>(lc.view().z().data());
      auto E = const_cast<float*>(lc.view().energy().data());
      std::unordered_map<float, std::vector<int>> map;

      for (int i = 0; i < lc->metadata().size(); ++i) {
        map[z[i]].push_back(i);
      }

      const int32_t n = static_cast<int32_t>(lc->metadata().size());
      std::cout << "Event Number of LCs " << n << std::endl;
      for (const auto& [Z, indices] : map) {
        std::cout << "z = " << Z << " -> Clusters : ";
        for (auto i : indices)
          std::cout << "\t( " << x[i] << ", " << y[i] << ", " << z[i] << ", " << E[i] << ")" << std::endl;
        std::cout << std::endl;
      }
      if (n > 0) {
        auto d_clIndex =
            cms::alpakatools::make_device_buffer<int[]>(queue, n);  // temporary buffer needed by CLUEstering
        auto dp_clIndex = const_cast<int*>(d_clIndex.data());
        clue::PointsDevice<3> d_points(queue, n, x, y, z, E, dp_clIndex);
        //  for(int iLC = 0; iLC < n; ++iLC){
        //  std::cout << "( " << x[iLC] << ", " << y[iLC] << ", " << z[iLC] << ", " << E[iLC] << " )" << std::endl;
        //  }
        //          auto isSeed =
        //            cms::alpakatools::make_device_buffer<int[]>(queue, nLCs);  // temporary buffer needed by CLUEstering

        clue::Clusterer<3> algo(queue, dc_, rho_, dm_);
        algo.make_clusters(queue, d_points);
        // create hosts points and do the copy
        clue::PointsHost<3> h_points(queue, n);
        clue::copyToHost(queue, h_points, d_points);
        alpaka::wait(queue);

        // get LCs indices in tracksters and fill the trackster collection
        const auto tsMap = clue::get_clusters(h_points);
        auto tracksters = std::vector<ticl::Trackster>(tsMap.size());
        for (long unsigned int i = 0; i < tsMap.size(); ++i) {
          const auto [beginLC, endLC] = tsMap.equal_range(i);
          std::copy(beginLC, endLC, std::back_inserter(tracksters[i].vertices()));
          tracksters[i].vertex_multiplicity().resize(tracksters[i].vertices().size(), 1);
        }

        std::cout << "Event Number of Tracksters " << tsMap.size() << std::endl;

        alpaka::memcpy(queue,
                       cms::alpakatools::make_host_view(h_points.view().coords[0], n),
                       cms::alpakatools::make_device_view(alpaka::getDev(queue), x, n),
                       (unsigned int)n);
        alpaka::memcpy(queue,
                       cms::alpakatools::make_host_view(h_points.view().coords[1], n),
                       cms::alpakatools::make_device_view(alpaka::getDev(queue), y, n),
                       (unsigned int)n);
        alpaka::memcpy(queue,
                       cms::alpakatools::make_host_view(h_points.view().coords[2], n),
                       cms::alpakatools::make_device_view(alpaka::getDev(queue), z, n),
                       (unsigned int)n);
        alpaka::memcpy(queue,
                       cms::alpakatools::make_host_view(h_points.weights(), n),
                       cms::alpakatools::make_device_view(alpaka::getDev(queue), E, n),
                       (unsigned int)n);
        alpaka::wait(queue);

        // compute trackster properties
        // TODO: merge with previous loop
        bool energyWeight = true;
        auto xHost = h_points.coords(0).data();
        auto yHost = h_points.coords(1).data();
        auto zHost = h_points.coords(2).data();
        auto EHost = h_points.weights();
        for (auto& trackster : tracksters) {
          size_t N = trackster.vertices().size();
          if (N == 0)
            continue;

          Eigen::Vector3f point;
          point << 0., 0., 0.;
          Eigen::Vector3f barycenter;
          barycenter << 0., 0., 0.;

          auto fillPoint = [&](const float x, const float y, const float z, const float weight = 1.f) {
            point[0] = weight * x;
            point[1] = weight * y;
            point[2] = weight * z;
          };

          // Initialize this trackster with default, dummy values
          trackster.setRawEnergy(0.f);
          trackster.setRawEmEnergy(0.f);
          trackster.setRawPt(0.f);
          trackster.setRawEmPt(0.f);

          float weight = 1.f / N;

          std::vector<float> layerClusterEnergies;

          for (size_t i = 0; i < N; ++i) {
            auto lcIdx = trackster.vertices(i);
            auto fraction = 1.f / trackster.vertex_multiplicity(i);
            trackster.addToRawEnergy(EHost[lcIdx] * fraction);
            // trackster.addToRawEmEnergy(EHost[lcIdx] * fraction);

            // Compute the weighted barycenter.
            if (energyWeight)
              weight = EHost[lcIdx] * fraction;
            fillPoint(xHost[lcIdx], yHost[lcIdx], zHost[lcIdx], weight);
            for (size_t j = 0; j < 3; ++j)
              barycenter[j] += point[j];

            layerClusterEnergies.push_back(EHost[lcIdx]);
          }
          float raw_energy = trackster.raw_energy();
          float inv_raw_energy = 1.f / raw_energy;
          if (energyWeight)
            barycenter *= inv_raw_energy;
          trackster.setBarycenter(ticl::Trackster::Vector(barycenter));

          trackster.calculateRawPt();
          trackster.calculateRawEmPt();

          std::cout << "  LC in TS: ";
          for (const auto& lc : trackster.vertices())
            std::cout << lc << " ";
          std::cout << std::endl;
          std::cout << "  energy raw: " << trackster.raw_energy() << std::endl;
          std::cout << "  barycenter: " << trackster.barycenter().x() << ", " << trackster.barycenter().y() << ", "
                    << trackster.barycenter().z() << std::endl;
        }
        // TODO: do a kernel that computes raw energy and barycenter and returns them

        iEvent.emplace(legacyTrackstersToken_, std::move(tracksters));
      } else {
        auto tracksters = std::vector<ticl::Trackster>();
        iEvent.emplace(legacyTrackstersToken_, std::move(tracksters));
      }
    }

  private:
    device::EDGetToken<HGCalSoAClustersDeviceCollection> const deviceTokenSoAClusters_;
    edm::EDPutTokenT<std::vector<ticl::Trackster>> const legacyTrackstersToken_;
    float rho_;
    std::array<float, 3> dc_;
    std::array<float, 3> dm_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(HeterogeneousTracksterProducer);
