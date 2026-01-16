#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HGCalReco/interface/Common.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "RecoHGCal/TICL/interface/TracksterLinkingAlgoBase.h"
#include "RecoHGCal/TICL/plugins/TracksterLinkingbyLayerOverlap.h"
#include "TICLGraph.h"
#include <numeric>
#include <algorithm>

using namespace ticl;

TracksterLinkingbyLayerOverlap::TracksterLinkingbyLayerOverlap(const edm::ParameterSet& conf,
                                                               edm::ConsumesCollector iC,
                                                               cms::Ort::ONNXRuntime const* onnxRuntime)
    : TracksterLinkingAlgoBase(conf, iC),
      min_trackster_energy_(conf.getParameter<double>("min_trackster_energy")),
      max_search_window_dR_(conf.getParameter<double>("max_search_window_dR")),
      max_layer_gap_(conf.getParameter<std::vector<int>>("max_layer_gap")),
      max_barycenter_dR_(conf.getParameter<std::vector<double>>("max_barycenter_dR")),
      min_pca_alignment_(conf.getParameter<double>("min_pca_alignment")),
      min_pca_quality_(conf.getParameter<double>("min_pca_quality")),
      max_sigma_timing_(conf.getParameter<double>("max_sigma_timing")),
      max_time_error_(conf.getParameter<double>("max_time_error")),
      require_timing_for_loose_geo_(conf.getParameter<bool>("require_timing_for_loose_geo")),
      w_layer_(conf.getParameter<double>("w_layer")),
      w_geo_(conf.getParameter<double>("w_geo")),
      w_pca_(conf.getParameter<double>("w_pca")),
      w_time_(conf.getParameter<double>("w_time")),
      stage1_min_energy_(conf.getParameter<double>("stage1_min_energy")),
      stage2_max_energy_(conf.getParameter<double>("stage2_max_energy")),
      stage2_min_pca_(conf.getParameter<double>("stage2_min_pca")) {}

void TracksterLinkingbyLayerOverlap::initialize(const HGCalDDDConstants* hgcons,
                                                const hgcal::RecHitTools rhtools,
                                                const edm::ESHandle<MagneticField> bfieldH,
                                                const edm::ESHandle<Propagator> propH) {
  hgcons_ = hgcons;
  rhtools_ = rhtools;
  bfield_ = bfieldH;
  propagator_ = propH;

  // Cache EM-Had interface position
  lastLayerEE_z_ = std::abs(rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z());
}

TracksterLinkingbyLayerOverlap::TracksterFeatures TracksterLinkingbyLayerOverlap::extractFeatures(
    const Trackster& ts, const std::vector<reco::CaloCluster>& layerClusters) {
  TracksterFeatures feat;

  // Basic properties
  feat.energy = ts.raw_energy();
  feat.barycenter = ts.barycenter();
  feat.eta = ts.barycenter().eta();
  feat.phi = ts.barycenter().phi();
  feat.pca_axis = ts.eigenvectors(0);  // Primary eigenvector

  // PCA quality
  auto const& eigenvalues = ts.eigenvalues();
  float sum = std::accumulate(std::begin(eigenvalues), std::end(eigenvalues), 0.f);
  feat.pca_quality = (sum > 0.f) ? eigenvalues[0] / sum : 0.f;

  // Layer extent
  feat.min_layer = std::numeric_limits<int>::max();
  feat.max_layer = std::numeric_limits<int>::min();

  for (auto lc_idx : ts.vertices()) {
    int layer = rhtools_.getLayerWithOffset(layerClusters[lc_idx].hitsAndFractions()[0].first);
    feat.min_layer = std::min(feat.min_layer, layer);
    feat.max_layer = std::max(feat.max_layer, layer);
  }

  // Timing (from trackster)
  feat.time = ts.time();
  feat.timeError = ts.timeError();
  feat.has_valid_timing = (feat.timeError > 0.f && feat.time != -99.f);

  return feat;
}

bool TracksterLinkingbyLayerOverlap::checkTimingCompatibility(const TracksterFeatures& inner,
                                                              const TracksterFeatures& outer,
                                                              LinkCandidate& candidate) {
  // If either trackster has invalid timing, can't check
  if (!inner.has_valid_timing || !outer.has_valid_timing) {
    candidate.timing_compatible = false;
    candidate.delta_time = -99.f;
    return false;
  }

  // Compute timing difference in units of sigma
  float deltaT = std::abs(inner.time - outer.time);
  float sigmaT = std::sqrt(inner.timeError * inner.timeError + outer.timeError * outer.timeError);

  // Guard against extremely uncertain timing (>max_time_error_ ns error)
  // CRITICAL FIX: If timing uncertainty is too large, it's meaningless for PU rejection
  if (sigmaT > max_time_error_) {
    candidate.timing_compatible = false;
    candidate.delta_time = -99.f;
    return false;
  }

  candidate.delta_time = deltaT / sigmaT;  // Store normalized delta

  // 3-sigma compatibility
  if (candidate.delta_time < max_sigma_timing_) {
    candidate.timing_compatible = true;
    return true;
  } else {
    candidate.timing_compatible = false;
    return false;
  }
}

bool TracksterLinkingbyLayerOverlap::checkCompatibility(const TracksterFeatures& feat_i,
                                                        const TracksterFeatures& feat_j,
                                                        LinkCandidate& candidate) {
  // Determine inner/outer by z-extent (choose lower min_layer as inner)
  const auto& ts_inner = (feat_i.min_layer <= feat_j.min_layer) ? feat_i : feat_j;
  const auto& ts_outer = (feat_i.min_layer <= feat_j.min_layer) ? feat_j : feat_i;

  candidate.inner_idx = ts_inner.trackster_idx;
  candidate.outer_idx = ts_outer.trackster_idx;

  // Check 1: Layer overlap/continuity
  int gap = std::abs(ts_outer.min_layer - ts_inner.max_layer);

  // Determine detector region for outer trackster (use appropriate layer gap threshold)
  bool is_outer_in_FH = (std::abs(ts_outer.barycenter.z()) > lastLayerEE_z_ + 50.0);  // FH starts ~50cm beyond CEH
  int max_gap = max_layer_gap_[is_outer_in_FH ? 1 : 0];

  if (gap > max_gap) {
    return false;  // Too far apart in layers
  }
  candidate.layer_gap = gap;

  // Check 2: Spatial window (coarse eta-phi cut)
  float dEta = std::abs(ts_inner.eta - ts_outer.eta);
  float dPhi = reco::deltaPhi(ts_inner.phi, ts_outer.phi);
  float dR_etaphi = std::sqrt(dEta * dEta + dPhi * dPhi);

  if (dR_etaphi > max_search_window_dR_) {
    return false;  // Outside search window
  }

  // Check 3: Barycenter alignment (projective distance)
  float dR_bary = reco::deltaR(
      ts_inner.barycenter.eta(), ts_inner.barycenter.phi(), ts_outer.barycenter.eta(), ts_outer.barycenter.phi());

  // Energy-dependent threshold (EE vs HAD)
  bool is_outer_in_HAD = (std::abs(ts_outer.barycenter.z()) > lastLayerEE_z_);
  float max_dR = max_barycenter_dR_[is_outer_in_HAD ? 1 : 0];

  if (dR_bary > max_dR) {
    return false;  // Barycenters too far apart
  }
  candidate.dR_bary = dR_bary;

  // Check 4: PCA alignment (NEW - currently unused in Skeletons!)
  float pca_dot = ts_inner.pca_axis.Dot(ts_outer.pca_axis);
  if (pca_dot < min_pca_alignment_) {
    return false;  // PCA axes not aligned (reject wide-angle merges)
  }
  candidate.pca_dot = pca_dot;

  // Check 5: Timing compatibility (NEW - critical for PU200!)
  if (!checkTimingCompatibility(ts_inner, ts_outer, candidate)) {
    // If timing fails and geometry is not very tight, reject
    if (require_timing_for_loose_geo_ && (dR_bary > 0.5 * max_dR || pca_dot < 0.95)) {
      return false;  // Loose geometry + bad timing → likely PU
    }
    // Otherwise allow (timing may be unreliable)
  }

  return true;  // All checks passed
}

float TracksterLinkingbyLayerOverlap::computeLinkScore(const LinkCandidate& candidate,
                                                       const TracksterFeatures& inner,
                                                       const TracksterFeatures& outer) {
  // Determine max gap for normalization
  bool is_outer_in_FH = (std::abs(outer.barycenter.z()) > lastLayerEE_z_ + 50.0);
  int max_gap = max_layer_gap_[is_outer_in_FH ? 1 : 0];

  // Layer continuity score (1.0 if gap=0, decreases with gap)
  float layer_score = 1.0 - (candidate.layer_gap / (float)max_gap);

  // Geometric proximity score (exponential decay with dR)
  float geo_score = std::exp(-candidate.dR_bary * candidate.dR_bary / 0.05);

  // PCA alignment score (linear with dot product)
  // CRITICAL FIX: Weight by minimum PCA quality to handle energy-asymmetric splits
  float pca_score = (candidate.pca_dot - min_pca_alignment_) / (1.0 - min_pca_alignment_);
  float min_pca_quality = std::min(inner.pca_quality, outer.pca_quality);
  pca_score *= min_pca_quality;  // Down-weight if either trackster has poor PCA

  // Timing score (1.0 if compatible, 0.0 if not, exponential decay)
  float time_score = 0.0;
  if (candidate.timing_compatible && candidate.delta_time >= 0.f) {
    time_score = std::exp(-candidate.delta_time * candidate.delta_time / 9.0);  // 9.0 = 3^2
  }

  // Weighted sum
  // CRITICAL FIX: Higher timing weight (0.3) at PU200 for better PU rejection
  float total_score = w_layer_ * layer_score + w_geo_ * geo_score + w_pca_ * pca_score + w_time_ * time_score;

  return total_score;
}

void TracksterLinkingbyLayerOverlap::buildLinkingGraph(std::vector<LinkCandidate>& candidates,
                                                       const std::vector<TracksterFeatures>& features,
                                                       std::vector<ticl::Node>& nodes) {
  // CRITICAL FIX: Deterministic sorting with tie-breaking
  // Use stable_sort + floating-point tolerance for determinism
  std::stable_sort(candidates.begin(), candidates.end(), [](const LinkCandidate& a, const LinkCandidate& b) {
    // Compare scores with floating-point tolerance
    if (std::abs(a.score - b.score) < 1e-6f) {
      // Deterministic tie-breaking: use trackster indices
      if (a.inner_idx != b.inner_idx) {
        return a.inner_idx < b.inner_idx;
      }
      return a.outer_idx < b.outer_idx;
    }
    return a.score > b.score;  // Higher score first
  });

  // Track which tracksters already have outgoing links
  std::vector<bool> has_outgoing(nodes.size(), false);

  // Iterate through candidates in score order
  for (const auto& cand : candidates) {
    unsigned int inner = cand.inner_idx;
    unsigned int outer = cand.outer_idx;

    // Only allow ONE outgoing link per trackster (avoid ambiguity)
    // But trackster can have MULTIPLE incoming links (merge many fragments)
    if (has_outgoing[inner]) {
      continue;  // Inner already linked to something better
    }

    // Create link: inner → outer
    nodes[inner].addOuterNeighbour(outer);
    nodes[outer].addInnerNeighbour(inner);

    has_outgoing[inner] = true;

    LogDebug("TracksterLinkingLayerOverlap")
        << "Linked trackster " << inner << " -> " << outer << " | score=" << cand.score
        << " | layer_gap=" << cand.layer_gap << " | dR=" << cand.dR_bary << " | PCA_dot=" << cand.pca_dot
        << " | timing=" << (cand.timing_compatible ? "OK" : "FAIL");
  }
}

void TracksterLinkingbyLayerOverlap::linkTracksters(
    const Inputs& input,
    std::vector<Trackster>& resultTracksters,
    std::vector<std::vector<unsigned int>>& linkedResultTracksters,
    std::vector<std::vector<unsigned int>>& linkedTracksterIdToInputTracksterId) {
  const auto& tracksters = input.tracksters;
  const auto& layerClusters = input.layerClusters;

  LogDebug("TracksterLinkingLayerOverlap")
      << "Starting Layer-Overlap linking for " << tracksters.size() << " tracksters";

  // Step 1: Extract features from all tracksters
  std::vector<TracksterFeatures> features;
  features.reserve(tracksters.size());

  for (size_t i = 0; i < tracksters.size(); ++i) {
    auto feat = extractFeatures(tracksters[i], layerClusters);
    feat.trackster_idx = i;

    // Filter by minimum energy and PCA quality
    if (feat.energy >= min_trackster_energy_ && feat.pca_quality >= min_pca_quality_) {
      features.push_back(feat);
    } else {
      LogDebug("TracksterLinkingLayerOverlap")
          << "Trackster " << i << " rejected: energy=" << feat.energy << " GeV, pca_quality=" << feat.pca_quality;
    }
  }

  LogDebug("TracksterLinkingLayerOverlap") << "After filtering: " << features.size() << " tracksters";

  // Step 2: Two-stage linking
  std::vector<LinkCandidate> stage1_candidates;  // High-confidence links
  std::vector<LinkCandidate> stage2_candidates;  // Satellite absorption

  // Stage 1: Link main fragments (both > stage1_min_energy)
  for (size_t i = 0; i < features.size(); ++i) {
    if (features[i].energy < stage1_min_energy_)
      continue;

    for (size_t j = i + 1; j < features.size(); ++j) {
      if (features[j].energy < stage1_min_energy_)
        continue;

      // Pre-filter by spatial proximity before expensive compatibility checks
      float dR = reco::deltaR(features[i].eta, features[i].phi, features[j].eta, features[j].phi);
      if (dR > max_search_window_dR_)
        continue;

      LinkCandidate candidate;
      if (checkCompatibility(features[i], features[j], candidate)) {
        candidate.score = computeLinkScore(candidate, features[candidate.inner_idx], features[candidate.outer_idx]);
        stage1_candidates.push_back(candidate);
      }
    }
  }

  LogDebug("TracksterLinkingLayerOverlap") << "Stage 1: " << stage1_candidates.size() << " candidate links";

  // Stage 2: Absorb small satellites
  for (size_t i = 0; i < features.size(); ++i) {
    if (features[i].energy < stage1_min_energy_)
      continue;  // Only large fragments can absorb

    for (size_t j = 0; j < features.size(); ++j) {
      if (i == j)
        continue;
      if (features[j].energy > stage2_max_energy_)
        continue;  // Only small fragments

      // Pre-filter by spatial proximity
      float dR = reco::deltaR(features[i].eta, features[i].phi, features[j].eta, features[j].phi);
      if (dR > max_search_window_dR_)
        continue;

      LinkCandidate candidate;
      if (checkCompatibility(features[i], features[j], candidate)) {
        // Stricter PCA requirement for satellites
        if (candidate.pca_dot >= stage2_min_pca_) {
          candidate.score = computeLinkScore(candidate, features[candidate.inner_idx], features[candidate.outer_idx]);
          candidate.score *= 0.8;  // Lower score for satellites to prefer stage 1 links
          stage2_candidates.push_back(candidate);
        }
      }
    }
  }

  LogDebug("TracksterLinkingLayerOverlap") << "Stage 2: " << stage2_candidates.size() << " satellite links";

  // Step 3: Combine candidates
  std::vector<LinkCandidate> all_candidates;
  all_candidates.reserve(stage1_candidates.size() + stage2_candidates.size());
  all_candidates.insert(all_candidates.end(), stage1_candidates.begin(), stage1_candidates.end());
  all_candidates.insert(all_candidates.end(), stage2_candidates.begin(), stage2_candidates.end());

  // Step 4: Build graph with best links only (avoid conflicts)
  std::vector<ticl::Node> nodes;
  nodes.reserve(tracksters.size());
  for (size_t i = 0; i < tracksters.size(); ++i) {
    nodes.emplace_back(i);
  }

  buildLinkingGraph(all_candidates, features, nodes);

  // Step 5: Extract connected components and merge
  TICLGraph graph(nodes);
  auto rootNodes = graph.getRootNodes();

  // Sort root nodes by energy for deterministic ordering
  std::sort(rootNodes.begin(), rootNodes.end(), [&tracksters](const ticl::Node& n1, const ticl::Node& n2) {
    unsigned int n1Id = n1.getId();
    unsigned int n2Id = n2.getId();
    return tracksters[n1Id].raw_energy() > tracksters[n2Id].raw_energy();
  });

  auto const& components = graph.findSubComponents(rootNodes);

  LogDebug("TracksterLinkingLayerOverlap") << "Found " << components.size() << " connected components";

  linkedTracksterIdToInputTracksterId.resize(components.size());

  int ic = 0;
  for (auto const& comp : components) {
    // Skip single-trackster components with very low energy
    if (comp.size() == 1) {
      if (tracksters[comp[0]].vertices().size() <= 3 && tracksters[comp[0]].raw_energy() < 5.f) {
        LogDebug("TracksterLinkingLayerOverlap") << "Skipping low-energy single trackster " << comp[0];
        continue;
      }
    }

    // Merge tracksters in this component
    Trackster outTrackster;
    std::vector<unsigned int> linkedTracksters;

    for (auto const& node : comp) {
      linkedTracksterIdToInputTracksterId[ic].push_back(node);
    }

    outTrackster.mergeTracksters(tracksters, linkedTracksterIdToInputTracksterId[ic]);
    linkedTracksters.push_back(resultTracksters.size());

    LogDebug("TracksterLinkingLayerOverlap") << "Component " << ic << ": merged " << comp.size()
                                             << " tracksters -> energy=" << outTrackster.raw_energy() << " GeV";

    resultTracksters.push_back(outTrackster);
    linkedResultTracksters.push_back(linkedTracksters);
    ++ic;
  }

  LogDebug("TracksterLinkingLayerOverlap") << "Output: " << resultTracksters.size() << " merged tracksters";
}
