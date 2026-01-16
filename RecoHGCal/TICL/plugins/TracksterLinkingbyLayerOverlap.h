#ifndef RecoHGCal_TICL_TracksterLinkingByLayerOverlap_H
#define RecoHGCal_TICL_TracksterLinkingByLayerOverlap_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "RecoHGCal/TICL/interface/TracksterLinkingAlgoBase.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TICLGraph.h"
#include <array>

namespace ticl {

  class TracksterLinkingbyLayerOverlap : public TracksterLinkingAlgoBase {
  public:
    TracksterLinkingbyLayerOverlap(const edm::ParameterSet& conf,
                                    edm::ConsumesCollector iC,
                                    cms::Ort::ONNXRuntime const* onnxRuntime = nullptr);

    ~TracksterLinkingbyLayerOverlap() override {}

    void linkTracksters(const Inputs& input,
                        std::vector<Trackster>& resultTracksters,
                        std::vector<std::vector<unsigned int>>& linkedResultTracksters,
                        std::vector<std::vector<unsigned int>>& linkedTracksterIdToInputTracksterId) override;

    void initialize(const HGCalDDDConstants* hgcons,
                    const hgcal::RecHitTools rhtools,
                    const edm::ESHandle<MagneticField> bfieldH,
                    const edm::ESHandle<Propagator> propH) override;

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc) {
      // Quality thresholds
      iDesc.add<double>("min_trackster_energy", 5.0);
      iDesc.add<double>("min_pca_quality", 0.70);

      // Spatial search
      iDesc.add<double>("max_search_window_dR", 0.15);

      // Layer overlap [EE+CEH, FH]
      iDesc.add<std::vector<int>>("max_layer_gap", {5, 8});

      // Geometric compatibility [EE, HAD]
      iDesc.add<std::vector<double>>("max_barycenter_dR", {0.08, 0.12});
      iDesc.add<double>("min_pca_alignment", 0.85);

      // Timing gates
      iDesc.add<double>("max_sigma_timing", 3.0);
      iDesc.add<double>("max_time_error", 10.0);
      iDesc.add<bool>("require_timing_for_loose_geo", true);

      // Scoring weights
      iDesc.add<double>("w_layer", 0.25);
      iDesc.add<double>("w_geo", 0.25);
      iDesc.add<double>("w_pca", 0.20);
      iDesc.add<double>("w_time", 0.30);

      // Multi-stage linking
      iDesc.add<double>("stage1_min_energy", 15.0);
      iDesc.add<double>("stage2_max_energy", 10.0);
      iDesc.add<double>("stage2_min_pca", 0.90);

      TracksterLinkingAlgoBase::fillPSetDescription(iDesc);
    }

    // Helper structs (public for use in implementation)
    struct TracksterFeatures {
      unsigned int trackster_idx;
      int min_layer;
      int max_layer;
      float energy;
      float eta;
      float phi;
      ticl::Vector barycenter;
      ticl::Vector pca_axis;
      float pca_quality;
      float time;
      float timeError;
      bool has_valid_timing;
    };

    struct LinkCandidate {
      unsigned int inner_idx;
      unsigned int outer_idx;
      float score;
      float layer_gap;
      float dR_bary;
      float pca_dot;
      float delta_time;
      bool timing_compatible;
    };

  private:
    using Vector = ticl::Trackster::Vector;

    // Core algorithm functions
    TracksterFeatures extractFeatures(const Trackster& ts,
                                       const std::vector<reco::CaloCluster>& layerClusters);

    bool checkCompatibility(const TracksterFeatures& feat_i,
                            const TracksterFeatures& feat_j,
                            LinkCandidate& candidate);

    bool checkTimingCompatibility(const TracksterFeatures& inner,
                                   const TracksterFeatures& outer,
                                   LinkCandidate& candidate);

    float computeLinkScore(const LinkCandidate& candidate,
                           const TracksterFeatures& inner,
                           const TracksterFeatures& outer);

    void buildLinkingGraph(std::vector<LinkCandidate>& candidates,
                           const std::vector<TracksterFeatures>& features,
                           std::vector<ticl::Node>& nodes);

    // Parameters
    float min_trackster_energy_;
    float max_search_window_dR_;
    std::vector<int> max_layer_gap_;  // [EE+CEH, FH]
    std::vector<double> max_barycenter_dR_;  // [EE, HAD]
    float min_pca_alignment_;
    float min_pca_quality_;
    float max_sigma_timing_;
    float max_time_error_;
    bool require_timing_for_loose_geo_;

    // Scoring weights
    float w_layer_;
    float w_geo_;
    float w_pca_;
    float w_time_;

    // Multi-stage parameters
    float stage1_min_energy_;
    float stage2_max_energy_;
    float stage2_min_pca_;

    // Geometry tools
    const HGCalDDDConstants* hgcons_;
    hgcal::RecHitTools rhtools_;
    float lastLayerEE_z_;  // Cached EM-Had interface position

    edm::ESHandle<MagneticField> bfield_;
    edm::ESHandle<Propagator> propagator_;
  };

}  // namespace ticl

#endif
