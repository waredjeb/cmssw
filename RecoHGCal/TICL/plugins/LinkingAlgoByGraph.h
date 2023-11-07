#ifndef RecoHGCal_TICL_LinkingAlgoByGraph_H__
#define RecoHGCal_TICL_LinkingAlgoByGraph_H__

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
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

namespace ticl {
  class LinkingAlgoByGraph final : public LinkingAlgoBase {
  public:
    LinkingAlgoByGraph(const edm::ParameterSet &conf);
    ~LinkingAlgoByGraph() override;

    void initialize(const HGCalDDDConstants *hgcons,
                    const hgcal::RecHitTools rhtools,
                    const edm::ESHandle<MagneticField> bfieldH,
                    const edm::ESHandle<Propagator> propH) override;

    void linkTracksters(const edm::Handle<std::vector<reco::Track>>,
                        const edm::Handle<std::vector<reco::CaloCluster>>, 
                        const edm::ValueMap<float> &,
                        const edm::ValueMap<float> &,
                        const edm::ValueMap<float> &,
                        const std::vector<reco::Muon> &,
                        const edm::Handle<std::vector<Trackster>>,
                        std::vector<TICLCandidate> &,
                        std::vector<TICLCandidate> &) override;

    float findSkeletonPoints(float percentage,
                            const float trackster_energy,
                            const std::vector<unsigned int> vertices,
                            const hgcal::RecHitTools &rhtools,
                            const std::vector<reco::CaloCluster> &layerClusters); 

    static void fillPSetDescription(edm::ParameterSetDescription &desc);

  private:
    using Vector = ticl::Trackster::Vector;

    void buildLayers();

    void dumpLinksFound(std::vector<std::vector<unsigned>> &resultCollection, const char *label) const;

    const float timing_quality_threshold_;
    const float del_;
    const float angle_first_cone_;
    const float angle_second_cone_;
    const float angle_third_cone_;
    const float max_height_cone_;

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
