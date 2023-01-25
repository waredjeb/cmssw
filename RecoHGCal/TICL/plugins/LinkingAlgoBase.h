#ifndef RecoHGCal_TICL_LinkingAlgoBase_H__
#define RecoHGCal_TICL_LinkingAlgoBase_H__

#include <vector>
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/HGCalReco/interface/EnergyRegressionAndIDModel.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "DataFormats/Math/interface/Vector3D.h"

namespace edm {
  class Event;
  class EventSetup;
}  // namespace edm

namespace ticl {
  class LinkingAlgoBase {
  public:
    LinkingAlgoBase(const edm::ParameterSet &conf) : algo_verbosity_(conf.getParameter<int>("algo_verbosity")) {}

    virtual ~LinkingAlgoBase(){};

    virtual void initialize(const HGCalDDDConstants *hgcons,
                            const hgcal::RecHitTools rhtools,
                            const edm::ESHandle<MagneticField> bfieldH,
                            const edm::ESHandle<Propagator> propH) = 0;

    virtual void linkTracksters(const edm::Handle<std::vector<reco::Track>>,
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
                                std::vector<bool>& masked_track) = 0;

    static void fillPSetDescription(edm::ParameterSetDescription &desc) { desc.add<int>("algo_verbosity", 0); };

    enum VerbosityLevel { None = 0, Basic, Advanced, Expert, Guru };

  protected:
    int algo_verbosity_;
  };
}  // namespace ticl

#endif
