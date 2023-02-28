
#include <cmath>
#include <string>

#include "RecoHGCal/TICL/plugins/SmoothingAlgoByMLP.h"

#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Common.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"

#include <bits/stdc++.h>

using namespace std;
using namespace ticl;
using namespace cms::Ort;


SmoothingAlgoByMLP::SmoothingAlgoByMLP(const edm::ParameterSet &conf)
    : LinkingAlgoBase(conf),
      del_tk_ts_layer1_(conf.getParameter<double>("delta_tk_ts_layer1")),
      del_tk_ts_int_(conf.getParameter<double>("delta_tk_ts_interface")),
      del_ts_em_had_(conf.getParameter<double>("delta_ts_em_had")),
      del_ts_had_had_(conf.getParameter<double>("delta_ts_had_had")),
      timing_quality_threshold_(conf.getParameter<double>("track_time_quality_threshold")),
      cutTk_(conf.getParameter<std::string>("cutTk")) {}

SmoothingAlgoByMLP::~SmoothingAlgoByMLP() {}

void SmoothingAlgoByMLP::initialize(
    const HGCalDDDConstants *hgcons,
    const hgcal::RecHitTools rhtools,
    const edm::ESHandle<MagneticField> bfieldH,
    const edm::ESHandle<Propagator> propH)
{
  hgcons_ = hgcons;
  rhtools_ = rhtools;

  bfield_ = bfieldH;
  propagator_ = propH;
}

void SmoothingAlgoByMLP::linkTracksters(
    const edm::Handle<std::vector<reco::Track>> tkH,
    const edm::ValueMap<float> &tkTime,
    const edm::ValueMap<float> &tkTimeErr,
    const edm::ValueMap<float> &tkTimeQual,
    const std::vector<reco::Muon> &muons,
    const edm::Handle<std::vector<Trackster>> tsH,
    std::vector<TICLCandidate> &resultLinked,
    std::vector<TICLCandidate> &chargedHadronsFromTk,
    std::vector<double>& prop_tracks_x,
    std::vector<double>& prop_tracks_y,
    std::vector<double>& prop_tracks_z,
    std::vector<double>& prop_tracks_eta,
    std::vector<double>& prop_tracks_phi,
    std::vector<double>& prop_tracks_px,
    std::vector<double>& prop_tracks_py,
    std::vector<double>& prop_tracks_pz,
    std::vector<bool>& masked_tracks,
    const TICLGraph &ticlGraph,
    const std::vector<reco::CaloCluster>& layerClusters,
    const ONNXRuntime *cache
) {
  const auto &tracksters = *tsH;
  long int N = tracksters.size();

  const float classification_threshold = 0; // no sigmoid, 0 = 0.5

  // CONFIGURATION OPTIONS
  const float radius = 30;
  const float energy_threshold = 10;

  /** PREPARING FEATURES **/

  const std::vector<std::string> input_names = {"features"};
  std::vector<float> features;

  const auto shapeFeatures = 63;

  std::vector<std::pair<unsigned, unsigned>> pairs;

  // Assuming this method is called per event
  // steps:
  // 2. get_trackster_representative_points (min z-point, max z-point)

  for (unsigned i = 0; i < tracksters.size(); ++i) {

    const auto &ts = tracksters[i];
    const float raw_energy = ts.raw_energy();

    // 1. we got a major trackster we want to smooth
    // ignore low energy tracksters
    if (raw_energy < energy_threshold) {
      continue;
    }

    const Vector &barycenter = ts.barycenter();
    const Vector &eigenvector0 = ts.eigenvectors(0);
    const std::array<float, 3> &eigenvalues = ts.eigenvalues();
    const std::array<float, 3> &sigmasPCA = ts.sigmasPCA();

    // 2. get representative points of the trackster (where (0, 0, 0) -> (bx, by, bz) intersects the min and max layer)
    const std::vector<unsigned int> &vertices_indices = ts.vertices();

    const auto max_z_lc_it = std::max_element(
      vertices_indices.begin(),
      vertices_indices.end(),
      [&layerClusters](const int &a, const int &b) {
        return layerClusters[a].z() > layerClusters[b].z();
      }
    );

    const auto min_z_lc_it = std::min_element(
      vertices_indices.begin(),
      vertices_indices.end(),
      [&layerClusters](const int &a, const int &b) {
        return layerClusters[a].z() > layerClusters[b].z();
      }
    );

    const reco::CaloCluster &min_z_lc = layerClusters[*min_z_lc_it];
    const reco::CaloCluster &max_z_lc = layerClusters[*max_z_lc_it];

    // compute the cylinder bounds
    const float t_min = min_z_lc.z() / barycenter.z();
    const float t_max = max_z_lc.z() / barycenter.z();

    const Vector x1 = Vector(
      t_min * barycenter.x(),
      t_min * barycenter.y(),
      min_z_lc.z()
    );

    const Vector x2 = Vector(
      t_max * barycenter.x(),
      t_max * barycenter.y(),
      min_z_lc.z()
    );

    // Loop over tracksters and see if they are in the cone
    for (unsigned ci = 0; ci < tracksters.size(); ++ci) {

      // no self loops
      if (ci == i) {
        continue;
      }

      // candidate trackster
      const auto &ct = tracksters[ci];

      const Vector &c_barycenter = ct.barycenter();
      const Vector &c_eigenvector0 = ct.eigenvectors(0);
      const std::array<float, 3> &c_eigenvalues = ct.eigenvalues();
      const std::array<float, 3> &c_sigmasPCA = ct.sigmasPCA();

      const std::vector<unsigned int> &c_vertices_indices = ct.vertices();

      const auto c_max_z_lc_it = std::max_element(
        c_vertices_indices.begin(),
        c_vertices_indices.end(),
        [&layerClusters](const int &a, const int &b) {
          return layerClusters[a].z() > layerClusters[b].z();
        }
      );

      const auto c_min_z_lc_it = std::min_element(
        c_vertices_indices.begin(),
        c_vertices_indices.end(),
        [&layerClusters](const int &a, const int &b) {
          return layerClusters[a].z() > layerClusters[b].z();
        }
      );

      const reco::CaloCluster &c_min_z_lc = layerClusters[*c_min_z_lc_it];
      const reco::CaloCluster &c_max_z_lc = layerClusters[*c_max_z_lc_it];

      // compute the distance and position in the cone
      if (c_barycenter.z() < x1.z() - radius || c_barycenter.z() > x2.z() + radius) {
        // not in z-cone bounds
        continue;
      }

      // d = np.linalg.norm(np.cross(x0 - x1, x0 - x2)) / np.linalg.norm(x2 - x1)
      const float distance = (c_barycenter - x1).Cross(c_barycenter - x2).R() / (x2 - x1).R();

      if (distance > radius) {
        // not in cone
        continue;
      }

      features.insert(features.end(), {
        (float) barycenter.x(),     // 0
        (float) barycenter.y(),     // 1
        (float) barycenter.z(),     // 2
        (float) barycenter.eta(),   // 3
        (float) barycenter.phi(), // 4
        (float) raw_energy,         // 5
        (float) ts.raw_em_energy(), // 6
        (float) eigenvalues[0],     // 7
        (float) eigenvalues[1],     // 8
        (float) eigenvalues[2],     // 9
        (float) eigenvector0.x(),   // 10
        (float) eigenvector0.y(),   // 11
        (float) eigenvector0.z(),   // 12
        (float) sigmasPCA[0],       // 13
        (float) sigmasPCA[1],       // 14
        (float) sigmasPCA[2],       // 15
        (float) c_barycenter.x(),   // 16
        (float) c_barycenter.y(),   // 17
        (float) c_barycenter.z(),   // 18
        (float) c_barycenter.eta(), // 19
        (float) c_barycenter.phi(), // 20
        (float) ct.raw_energy(),    // 21
        (float) ct.raw_em_energy(), // 22
        (float) c_eigenvalues[0],   // 23
        (float) c_eigenvalues[1],   // 24
        (float) c_eigenvalues[2],   // 25
        (float) c_eigenvector0.x(), // 26
        (float) c_eigenvector0.y(), // 27
        (float) c_eigenvector0.z(), // 28
        (float) c_sigmasPCA[0],     // 29
        (float) c_sigmasPCA[1],     // 30
        (float) c_sigmasPCA[2],     // 31
        (float) min_z_lc.x(),       // 32
        (float) min_z_lc.y(),       // 33
        (float) min_z_lc.z(),       // 34
        (float) max_z_lc.x(),       // 35
        (float) max_z_lc.y(),       // 36
        (float) max_z_lc.z(),       // 37
        (float) c_min_z_lc.x(),     // 38
        (float) c_min_z_lc.y(),     // 39
        (float) c_min_z_lc.z(),     // 40
        (float) c_max_z_lc.x(),     // 41
        (float) c_max_z_lc.y(),     // 42
        (float) c_max_z_lc.z(),     // 43
        (float) ts.id_probabilities(0), // 44
        (float) ts.id_probabilities(1), // 45
        (float) ts.id_probabilities(2), // 46
        (float) ts.id_probabilities(3), // 47
        (float) ts.id_probabilities(4), // 48
        (float) ts.id_probabilities(5), // 49
        (float) ts.id_probabilities(6), // 50
        (float) ts.id_probabilities(7), // 51
        (float) ct.id_probabilities(0), // 52
        (float) ct.id_probabilities(1), // 53
        (float) ct.id_probabilities(2), // 54
        (float) ct.id_probabilities(3), // 55
        (float) ct.id_probabilities(4), // 56
        (float) ct.id_probabilities(5), // 57
        (float) ct.id_probabilities(6), // 58
        (float) ct.id_probabilities(7), // 59
        (float) distance,               // 60
        (float) ts.vertices().size(),   // 61
        (float) ct.vertices().size()    // 62
      });

      // keep track of the samples
      pairs.push_back(std::make_pair(i, ci));
    }
  }

  std::cout << "PAIRS EXTRACTED: " << pairs.size() << std::endl;

  /** RUNNING THE NETWORK **/
  std::vector<float> edge_predictions;
  if (pairs.size() > 0) {

    // Prepare network input
    std::vector<std::vector<int64_t>> input_shapes;

    FloatArrays data;
    data.emplace_back(features);

    input_shapes.push_back({ (int64_t) pairs.size(), shapeFeatures });

    // get the network output
    // input_names: list of the names of the input nodes ("features")
    // input_values: list of input arrays for each input node. The order of `input_values` must match `input_names`.
    // input_shapes: list of `int64_t` arrays specifying the shape of each input node. Can leave empty if the model does not have dynamic axes.
    // output_names: names of the output nodes to get outputs from. Empty list means all output nodes.
    // batch_size: number of samples in the batch. Each array in `input_values` must have a shape layout of (batch_size, ...).
    edge_predictions = cache->run(input_names, data, input_shapes, {}, (int64_t) pairs.size())[0];
  }

  // interpret the results
  std::vector<TICLCandidate> connectedCandidates;

  for (int trackster_id = 0; trackster_id < N; ++trackster_id) {
    bool skip = false;

    TICLCandidate tracksterCandidate;
    tracksterCandidate.addTrackster(edm::Ptr<Trackster>(tsH, trackster_id));

    int yi = 0;
    for (auto &p : pairs) {
      const int pi = p.first;
      const int pj = p.second;
      const float score = edge_predictions[yi++];

      if (pj == trackster_id && score > classification_threshold) {
        // the trackster is connected to another trackster
        skip = true;
        break;
      }

      if (trackster_id == pi && score > classification_threshold) {
        // trackster is the main trackster
        // check if the score is > threshold
        std::cout << "MERGING TRACKSTERS: (" << pi << ", " << pj << ")" << std::endl;
        tracksterCandidate.addTrackster(edm::Ptr<Trackster>(tsH, pj));
      }
    }

    if (!skip) {
      connectedCandidates.push_back(tracksterCandidate);
    }
  }

  std::cout << "MLP Smoothing: " << N << " -> " << connectedCandidates.size() << std::endl;

  // The final candidates are passed to `resultLinked`
  resultLinked.insert(std::end(resultLinked), std::begin(connectedCandidates), std::end(connectedCandidates));
}  // linkTracksters


void SmoothingAlgoByMLP::fillPSetDescription(edm::ParameterSetDescription &desc) {
  desc.add<std::string>("cutTk",
                        "1.48 < abs(eta) < 3.0 && pt > 1. && quality(\"highPurity\") && "
                        "hitPattern().numberOfLostHits(\"MISSING_OUTER_HITS\") < 5");
  desc.add<double>("delta_tk_ts_layer1", 0.02);
  desc.add<double>("delta_tk_ts_interface", 0.03);
  desc.add<double>("delta_ts_em_had", 0.03);
  desc.add<double>("delta_ts_had_had", 0.03);
  desc.add<double>("track_time_quality_threshold", 0.5);
  LinkingAlgoBase::fillPSetDescription(desc);
}
