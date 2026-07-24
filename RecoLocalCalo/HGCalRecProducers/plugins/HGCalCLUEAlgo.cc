#include "RecoLocalCalo/HGCalRecProducers/plugins/HGCalCLUEAlgo.h"

// Geometry
#include "DataFormats/CaloRecHit/interface/CaloClusterHostCollection.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/TICL/interface/AssociationMap.h"
#include "DataFormats/TICL/interface/FillAssociator.h"
#include "DataFormats/TICL/interface/HitsAndFractionsHost.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

#include "CLUEstering/CLUEstering.hpp"

#include <algorithm>
#include <limits>
#include <numeric>
#include <span>

using namespace hgcal_clustering;

template <typename T, typename STRATEGY>
void HGCalCLUEAlgoT<T, STRATEGY>::getEventSetupPerAlgorithm(const edm::EventSetup& es) {
  cells_.clear();
  numberOfClustersPerLayer_.clear();
  cells_.resize(2 * (maxlayer_ + 1));
  numberOfClustersPerLayer_.resize(2 * (maxlayer_ + 1), 0);
}

template <typename T, typename STRATEGY>
void HGCalCLUEAlgoT<T, STRATEGY>::populate(const HGCRecHitCollection& hits) {
  // loop over all hits and create the Hexel structure, skip energies below ecut
  if (dependSensor_) {
    // for each layer and wafer calculate the thresholds (sigmaNoise and energy)
    // once
    computeThreshold();
  }

  for (unsigned int i = 0; i < hits.size(); ++i) {
    const HGCRecHit& hgrh = hits[i];
    DetId detid = hgrh.detid();
    unsigned int layerOnSide = (rhtools_.getLayerWithOffset(detid) - 1);

    // set sigmaNoise default value 1 to use kappa value directly in case of
    // sensor-independent thresholds
    float sigmaNoise = 1.f;
    if (dependSensor_) {
      int thickness_index = rhtools_.getSiThickIndex(detid);
      if (thickness_index == -1)
        thickness_index = maxNumberOfThickIndices_;

      double storedThreshold = thresholds_[layerOnSide][thickness_index];
      if (detid.det() == DetId::HGCalHSi || detid.subdetId() == HGCHEF) {
        storedThreshold = thresholds_[layerOnSide][thickness_index + deltasi_index_regemfac_];
      }
      sigmaNoise = v_sigmaNoise_[layerOnSide][thickness_index];

      if (hgrh.energy() < storedThreshold)
        continue;  // this sets the ZS threshold at ecut times the sigma noise
                   // for the sensor
    }
    if (!dependSensor_ && hgrh.energy() < ecut_)
      continue;
    const GlobalPoint position(rhtools_.getPosition(detid));
    int offset = ((rhtools_.zside(detid) + 1) >> 1) * maxlayer_;
    int layer = layerOnSide + offset;
    // setting the layer position only once per layer
    if (cells_[layer].layerDim3 == std::numeric_limits<float>::infinity())
      cells_[layer].layerDim3 = position.z();

    cells_[layer].detid.emplace_back(detid);
    if constexpr (std::is_same_v<STRATEGY, HGCalScintillatorStrategy>) {
      cells_[layer].dim1.emplace_back(position.eta());
      cells_[layer].dim2.emplace_back(position.phi());
    }  // else, isSilicon == true and eta phi values will not be used
    else {
      cells_[layer].dim1.emplace_back(position.x());
      cells_[layer].dim2.emplace_back(position.y());
    }
    cells_[layer].weight.emplace_back(hgrh.energy());
    cells_[layer].sigmaNoise.emplace_back(sigmaNoise);
  }
}

// Run the external CLUEstering (2D, per layer) algorithm to populate, for every
// layer, the per-hit cluster index and the list of seeds, together with the
// number of clusters found on the layer.
template <typename T, typename STRATEGY>
void HGCalCLUEAlgoT<T, STRATEGY>::makeClusters() {
  for (auto l = 0u; l < 2 * maxlayer_ + 2; ++l) {
    // The critical distance for the local-density calculation (vecDeltasC_), the
    // seed-promotion distance (vecDeltasSeed_) and the outlier distance
    // (vecDeltasO_) are three independent, per-sensor-group parameters.
    unsigned int deltaIdx = 3;  // scintillator sensor group
    if constexpr (std::is_same_v<STRATEGY, HGCalSiliconStrategy>) {
      if (l % maxlayer_ < lastLayerEE_)
        deltaIdx = 0;
      else if (l % maxlayer_ < (firstLayerBH_ - 1))
        deltaIdx = 1;
      else
        deltaIdx = 2;
    }
    const float delta_c = vecDeltasC_[deltaIdx];
    const float delta_seed = vecDeltasSeed_[deltaIdx];
    const float delta_o = vecDeltasO_[deltaIdx];

    const auto nhits = cells_[l].dim1.size();
    cells_[l].clusterIndex.resize(nhits, -1);
    if (nhits == 0u)
      continue;

    auto queue = clue::get_queue(0u);
    auto clusterer = clue::Clusterer<2>(delta_c, kappa_, delta_o, delta_seed);
    auto points = clue::PointsHost<2>(
        queue, static_cast<int32_t>(nhits), cells_[l].dim1, cells_[l].dim2, cells_[l].weight, cells_[l].clusterIndex);
    points.set_density_uncertainty(std::span<float>(cells_[l].sigmaNoise));
    clusterer.make_clusters(points);
    numberOfClustersPerLayer_[l] = points.n_clusters();

    const auto seeds = clusterer.getSeeds();
    std::ranges::copy(seeds, std::back_inserter(cells_[l].seeds));
  }
#if DEBUG_CLUSTERS_ALPAKA
  hgcalUtils::DumpLegacySoA dumperLegacySoA;
  dumperLegacySoA.dumpInfos(cells_, moduleType_);
#endif
}

template <typename T, typename STRATEGY>
std::vector<reco::BasicCluster> HGCalCLUEAlgoT<T, STRATEGY>::getClustersLegacy(bool) {
  // The legacy std::vector<reco::CaloCluster> product is built by the producer
  // from the SoA and the transient hits-and-fractions association map.
  return std::vector<reco::BasicCluster>();
}

template <typename T, typename STRATEGY>
ticl::LayerClustersAndAssociations HGCalCLUEAlgoT<T, STRATEGY>::getClusters(bool) {
  std::vector<int> offsets(numberOfClustersPerLayer_.size(), 0);
  for (unsigned layerId = 1; layerId < offsets.size(); ++layerId) {
    offsets[layerId] = offsets[layerId - 1] + numberOfClustersPerLayer_[layerId - 1];
  }
  const auto totalNumberOfClusters = offsets.back() + numberOfClustersPerLayer_.back();

  const auto total_rechits = std::accumulate(
      cells_.begin(), cells_.end(), 0, [](auto acc, const auto& cell) { return acc + cell.dim1.size(); });

  ticl::LayerClustersAndAssociations clusters_and_associations(totalNumberOfClusters, total_rechits);
  auto layer_clusters_view = clusters_and_associations.layer_clusters->view();

  // keys (per clustered hit -> global cluster index) and values ({DetId, fraction})
  // are filled in the same loop so that they stay aligned entry-by-entry.
  std::vector<int> cluster_hit_associations;
  std::vector<ticl::HitAndFraction> detid_and_fractions;

  for (unsigned int layerId = 0; layerId < 2 * maxlayer_ + 2; ++layerId) {
    if (cells_[layerId].dim1.empty() || numberOfClustersPerLayer_[layerId] == 0)
      continue;

    auto queue = clue::get_queue(0u);
    auto points = clue::make_clustered_points<2>(queue,
                                                 std::span<const float>(cells_[layerId].dim1),
                                                 std::span<const float>(cells_[layerId].dim2),
                                                 std::span<const float>(cells_[layerId].weight),
                                                 std::span<const int>(cells_[layerId].clusterIndex));
    if (points.size() <= 0)
      continue;

    const auto clusters = clue::get_clusters(points);
    const auto weights = points.weights();
    const auto coords0 = points.coords(0);
    const auto coords1 = points.coords(1);

    for (auto cl = 0u; cl < clusters.size(); ++cl) {
      const auto cluster = clusters[cl];  // span of point indices belonging to this cluster
      const auto globalClusterIndex = cl + offsets[layerId];

      float energy = 0.f;
      float maxEnergy = std::numeric_limits<float>::lowest();
      int maxEnergyIdx = cluster.front();
      for (auto p : cluster) {
        energy += weights[p];
        if (weights[p] > maxEnergy) {
          maxEnergy = weights[p];
          maxEnergyIdx = p;
        }
        cluster_hit_associations.push_back(globalClusterIndex);
        detid_and_fractions.push_back(ticl::HitAndFraction{cells_[layerId].detid[p], 1.f});
      }

      float x = 0.f;
      float y = 0.f;
      const float z = cells_[layerId].layerDim3;
      if constexpr (std::is_same_v<STRATEGY, HGCalSiliconStrategy>) {
        const auto max_energy_detid = cells_[layerId].detid[maxEnergyIdx];
        const auto thick = rhtools_.getSiThickIndex(max_energy_detid);
        float total_weight_log = 0.f;
        for (auto p : cluster) {
          const float d1 = coords0[p] - coords0[maxEnergyIdx];
          const float d2 = coords1[p] - coords1[maxEnergyIdx];
          if ((d1 * d1 + d2 * d2) < positionDeltaRho2_) {
            const float Wi = std::max(thresholdW0_[thick] + std::log(weights[p] / energy), 0.);
            x += coords0[p] * Wi;
            y += coords1[p] * Wi;
            total_weight_log += Wi;
          }
        }
        if (total_weight_log != 0.f) {
          const float inv_tot_weight = 1.f / total_weight_log;
          x *= inv_tot_weight;
          y *= inv_tot_weight;
        } else {
          x = coords0[maxEnergyIdx];
          y = coords1[maxEnergyIdx];
        }
      } else {
        const auto centroid = clue::weighted_cluster_centroid<2>(points, cl);
        x = centroid[0];
        y = centroid[1];
      }

      layer_clusters_view.position().x()[globalClusterIndex] = x;
      layer_clusters_view.position().y()[globalClusterIndex] = y;
      layer_clusters_view.position().z()[globalClusterIndex] = z;
      layer_clusters_view.position().layer()[globalClusterIndex] = static_cast<int>(layerId);
      layer_clusters_view.position().cells()[globalClusterIndex] = static_cast<int>(cluster.size());
      layer_clusters_view.energy().energy()[globalClusterIndex] = energy;
      layer_clusters_view.energy().correctedEnergy()[globalClusterIndex] = -1.f;
      layer_clusters_view.energy().correctedEnergyUncertainty()[globalClusterIndex] = -1.f;
      layer_clusters_view.indexes().caloID()[globalClusterIndex] = reco::CaloID(reco::CaloID::DET_HGCAL_ENDCAP);
      layer_clusters_view.indexes().algoID()[globalClusterIndex] = algoId_;
      layer_clusters_view.indexes().seedID()[globalClusterIndex] = cells_[layerId].detid[cells_[layerId].seeds[cl]];
      layer_clusters_view.indexes().flags()[globalClusterIndex] = 0;
    }
  }

  // Rebuild the transient hits-and-fractions association map with the exact
  // number of clustered hits and fill it with the (key, value) pairs above.
  auto new_hits_and_fractions = std::make_unique<ticl::HitsAndFractionsHost>(
      cms::alpakatools::host(), cluster_hit_associations.size(), totalNumberOfClusters);
  clusters_and_associations.hits_and_fractions = std::move(new_hits_and_fractions);

  alpaka_serial_sync::Queue queue(cms::alpakatools::host());
  ticl::associator::fill<alpaka_serial_sync::Acc1D>(
      queue,
      clusters_and_associations.hits_and_fractions->view(),
      static_cast<std::span<const int>>(cluster_hit_associations),
      static_cast<std::span<const ticl::HitAndFraction>>(detid_and_fractions));

  return clusters_and_associations;
}

template <typename T, typename STRATEGY>
void HGCalCLUEAlgoT<T, STRATEGY>::computeThreshold() {
  // To support the TDR geometry and also the post-TDR one (v9 onwards), we
  // need to change the logic of the vectors containing signal to noise and
  // thresholds. The first 3 indices will keep on addressing the different
  // thicknesses of the Silicon detectors in CE_E , the next 3 indices will
  // address the thicknesses of the Silicon detectors in CE_H, while the last
  // one, number 6 (the seventh) will address the Scintillators. This change
  // will support both geometries at the same time.

  if (initialized_)
    return;  // only need to calculate thresholds once

  initialized_ = true;

  std::vector<double> dummy;

  dummy.resize(maxNumberOfThickIndices_ + !isNose_,
               0);  // +1 to accomodate for the Scintillators
  thresholds_.resize(maxlayer_, dummy);
  v_sigmaNoise_.resize(maxlayer_, dummy);

  for (unsigned ilayer = 1; ilayer <= maxlayer_; ++ilayer) {
    for (unsigned ithick = 0; ithick < maxNumberOfThickIndices_; ++ithick) {
      float sigmaNoise = 0.001f * fcPerEle_ * nonAgedNoises_[ithick] * dEdXweights_[ilayer] /
                         (fcPerMip_[ithick] * thicknessCorrection_[ithick]);
      thresholds_[ilayer - 1][ithick] = sigmaNoise * ecut_;
      v_sigmaNoise_[ilayer - 1][ithick] = sigmaNoise;
      LogDebug("HGCalCLUEAlgo") << "ilayer: " << ilayer << " nonAgedNoises: " << nonAgedNoises_[ithick]
                                << " fcPerEle: " << fcPerEle_ << " fcPerMip: " << fcPerMip_[ithick]
                                << " noiseMip: " << fcPerEle_ * nonAgedNoises_[ithick] / fcPerMip_[ithick]
                                << " sigmaNoise: " << sigmaNoise << "\n";
    }

    if (!isNose_) {
      float scintillators_sigmaNoise = 0.001f * noiseMip_ * dEdXweights_[ilayer] / sciThicknessCorrection_;
      thresholds_[ilayer - 1][maxNumberOfThickIndices_] = ecut_ * scintillators_sigmaNoise;
      v_sigmaNoise_[ilayer - 1][maxNumberOfThickIndices_] = scintillators_sigmaNoise;
      LogDebug("HGCalCLUEAlgo") << "ilayer: " << ilayer << " noiseMip: " << noiseMip_
                                << " scintillators_sigmaNoise: " << scintillators_sigmaNoise << "\n";
    }
  }
}

// explicit template instantiation
template class HGCalCLUEAlgoT<HGCalSiliconLayerTiles, HGCalSiliconStrategy>;
template class HGCalCLUEAlgoT<HGCalScintillatorLayerTiles, HGCalScintillatorStrategy>;
template class HGCalCLUEAlgoT<HFNoseLayerTiles, HGCalSiliconStrategy>;
