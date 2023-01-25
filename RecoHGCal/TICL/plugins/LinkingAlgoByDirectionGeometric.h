#ifndef RecoHGCal_TICL_LinkingAlgoByDirectionGeometric_H__
#define RecoHGCal_TICL_LinkingAlgoByDirectionGeometric_H__

#include <memory>
#include <array>
#include "RecoHGCal/TICL/plugins/LinkingAlgoBase.h"
#include "RecoHGCal/TICL/interface/commons.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"

#include "PhysicsTools/TensorFlow/interface/TfGraphRecord.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "PhysicsTools/TensorFlow/interface/TfGraphDefWrapper.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

namespace ticl {
  class LinkingAlgoByDirectionGeometric final : public LinkingAlgoBase {
  public:
    LinkingAlgoByDirectionGeometric(const edm::ParameterSet &conf);
    ~LinkingAlgoByDirectionGeometric() override;

    void initialize(const HGCalDDDConstants *hgcons,
                    const hgcal::RecHitTools rhtools,
                    const edm::ESHandle<MagneticField> bfieldH,
                    const edm::ESHandle<Propagator> propH) override;

    void linkTracksters(const edm::Handle<std::vector<reco::Track>>,
                        const edm::ValueMap<float> &,
                        const edm::ValueMap<float> &,
                        const edm::ValueMap<float> &,
                        const std::vector<reco::Muon> &,
                        const edm::Handle<std::vector<Trackster>>,
                        const std::vector<reco::CaloCluster> &layerClusters,
                        const edm::ValueMap<std::pair<float, float>> &layerClustersTime,
                        std::vector<Trackster> &tracksterMergeCollectionResult,
                        std::vector<TICLCandidate> &,
                        std::vector<TICLCandidate> &,
                        const EnergyRegressionAndIDModel &,
                        std::vector<float> &,
                        std::vector<float> &,
                        std::vector<float> &,
                        std::vector<int> &,
                        std::vector<double>& prop_tracks_x,
                        std::vector<double>& prop_tracks_y,
                        std::vector<double>& prop_tracks_z,
                        std::vector<double>& prop_tracks_eta,
                        std::vector<double>& prop_tracks_phi,
                        std::vector<double>& prop_tracks_px,
                        std::vector<double>& prop_tracks_py,
                        std::vector<double>& prop_tracks_pz,
                        std::vector<bool>& masked_track) override;

      

    static void fillPSetDescription(edm::ParameterSetDescription &desc);

  private:
    typedef math::XYZVector Vector;
    typedef std::vector<double> Vec;

    void buildLayers();

    void fillTrackstersTile(const std::vector<Trackster> &t, std::array<TICLLayerTile, 2> &tracksterTiles);

    void fillTracksTile(const std::vector<std::pair<Vector, unsigned>> &propTracks,
                        std::array<TICLLayerTile, 2> &tracksTiles);

    void findSmallTrackstersInWindow(const std::vector<Trackster> &tracksters,
                                     const std::vector<int> &maskSmall,
                                     const std::array<TICLLayerTile, 2> &tracksterTiles,
                                     const float delta,
                                     const float separation,
                                     std::vector<std::vector<unsigned>> &resultCollection,
                                     std::vector<float> &distancesVec,
                                     std::vector<int> &distancesVecIdx,
                                     bool useMask);
    void findTrackstersInWindow(const std::vector<Trackster> &tracksters,
                                const std::array<TICLLayerTile, 2> &tracksterTiles,
                                const float delta,
                                const float separation,
                                std::vector<std::vector<unsigned>> &resultCollection,
                                std::vector<float> &distancesVec,
                                std::vector<int> &distancesVecIdx,
                                bool useMask);

    void tracksterToTrackLinking(std::vector<Trackster> &tracksterMergeCollection,
                                 std::vector<std::vector<unsigned>> &tracksterMergeCollectionIndices,
                                 const std::vector<Trackster> &trackster,
                                 const edm::Handle<std::vector<Trackster>> tsH,
                                 const std::vector<reco::Track> &tracks,
                                 const edm::Handle<std::vector<reco::Track>> tkH,
                                 const edm::ValueMap<float> &tkTime,
                                 const edm::ValueMap<float> &tkTimeErr,
                                 const edm::ValueMap<float> &tkTimeQual,
                                 std::vector<std::pair<Vector, unsigned>> &propTracks,
                                 std::array<TICLLayerTile, 2> &tracksterMergeTiles,
                                 const float delta,
                                 const float separation,
                                 std::vector<TICLCandidate> &candidate,
                                 std::vector<TICLCandidate> &chargedCandidatesFromTracks,
                                 std::vector<float> &,
                                 std::vector<float> &);

    bool timeCompatible(const Trackster &, const Trackster &);

    bool timeAndEnergyCompatible(const reco::Track &track,
                                 const Trackster &trackster,
                                 const float &tkTime,
                                 const float &tkTimeErr,
                                 const float &tkTimeQual);

    void recordTrackster(const unsigned ts,  // trackster index
                         const std::vector<Trackster> &tracksters,
                         const edm::Handle<std::vector<Trackster>> tsH,
                         std::vector<unsigned> &ts_mask,
                         float &energy_in_candidate,
                         TICLCandidate &candidate);

    void energyRegressionAndID(const std::vector<reco::CaloCluster> &layerClusters,
                               const tensorflow::Session *eidSession,
                               std::vector<Trackster> &tracksters) const;
    void dumpLinksFound(std::vector<std::vector<unsigned>> &resultCollection, const char *label) const;

    const float tkEnergyCut_ = 2.0f;
    const float maxDeltaT_ = 3.0f;
    const float del_tk_ts_layer1_;
    const float del_tk_ts_int_;
    const float del_ts_em_had_;
    const float del_ts_had_had_;
    const float separationSmall_threshold_;
    const float separation_threshold_;
    const int maxDepth_;
    const float timing_quality_threshold_;

    const StringCutObjectSelector<reco::Track> cutTk_;
    std::once_flag initializeGeometry_;

    const HGCalDDDConstants *hgcons_;

    std::unique_ptr<GeomDet> firstDisk_[2];
    std::unique_ptr<GeomDet> interfaceDisk_[2];

    hgcal::RecHitTools rhtools_;

    edm::ESHandle<MagneticField> bfield_;
    edm::ESHandle<Propagator> propagator_;
  };
}  // namespace ticl
#endif
