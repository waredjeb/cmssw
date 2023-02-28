
#ifndef RecoHGCal_TICL_LinkingAlgoByGNN_H__
#define RecoHGCal_TICL_LinkingAlgoByGNN_H__

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

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

using namespace cms::Ort;
namespace ticl {
  class LinkingAlgoByGNN final : public LinkingAlgoBase {
  public:
    LinkingAlgoByGNN(const edm::ParameterSet &conf);
    ~LinkingAlgoByGNN() override;

    void initialize(const HGCalDDDConstants *hgcons,
                    const hgcal::RecHitTools rhtools,
                    const edm::ESHandle<MagneticField> bfieldH,
                    const edm::ESHandle<Propagator> propH) override;

    void linkTracksters(const edm::Handle<std::vector<reco::Track>> tkH,
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
                                                     const TICLGraph &,
                                                     const std::vector<reco::CaloCluster> &,
                                                     const ONNXRuntime *cache) override;

    static void fillPSetDescription(edm::ParameterSetDescription &desc);

  private:
    typedef math::XYZVector Vector;

    void dumpLinksFound(std::vector<std::vector<unsigned>> &resultCollection, const char *label) const;

    const float tkEnergyCut_ = 2.0f;
    const float maxDeltaT_ = 3.0f;
    const float del_tk_ts_layer1_;
    const float del_tk_ts_int_;
    const float del_ts_em_had_;
    const float del_ts_had_had_;

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
