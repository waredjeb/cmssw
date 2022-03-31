#ifndef SimDataFormats_Associations_TracksterToSimTracksterHitLCAssociatorBaseImpl_h
#define SimDataFormats_Associations_TracksterToSimTracksterHitLCAssociatorBaseImpl_h

/** \class TracksterToSimTracksterHitLCAssociatorBaseImpl
 *
 * Base class for TracksterToSimTracksterHitLCAssociator. Methods take as input
 * the handles of Tracksters, LayerClusters and SimTracksters collections and return an
 * AssociationMap (oneToManyWithQuality)
 *
 *  \author Leonardo Cristella
 */

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

//#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
typedef std::vector<SimCluster> SimClusterCollection;
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"

namespace hgcal {

  enum validationType { Linking = 0, PatternRecognition, PatternRecognition_CP };

  typedef edm::AssociationMap<
      edm::OneToManyWithQualityGeneric<ticl::TracksterCollection, ticl::TracksterCollection, std::pair<float, float>>>
      SimToRecoCollectionSimTracksters;
  typedef edm::AssociationMap<
      edm::OneToManyWithQualityGeneric<ticl::TracksterCollection, ticl::TracksterCollection, float>>
      RecoToSimCollectionSimTracksters;

  class TracksterToSimTracksterHitLCAssociatorBaseImpl {
  public:
    /// Constructor
    TracksterToSimTracksterHitLCAssociatorBaseImpl();
    /// Destructor
    virtual ~TracksterToSimTracksterHitLCAssociatorBaseImpl();

    /// Associate a Trackster to SimClusters
    virtual hgcal::RecoToSimCollectionSimTracksters associateRecoToSim(
        const edm::Handle<ticl::TracksterCollection> &tCH,
        const edm::Handle<reco::CaloClusterCollection> &lCCH,
        const edm::Handle<SimClusterCollection> &sCCH,
        const edm::Handle<CaloParticleCollection> &cPCH,
        const edm::Handle<std::map<uint, std::vector<uint>>> &simTrackstersMapH,
        const edm::Handle<ticl::TracksterCollection> &sTCH,
        const edm::Handle<ticl::TracksterCollection> &sTfromCPCH,
        const validationType valType) const;

    /// Associate a SimCluster to Tracksters
    virtual hgcal::SimToRecoCollectionSimTracksters associateSimToReco(
        const edm::Handle<ticl::TracksterCollection> &tCH,
        const edm::Handle<reco::CaloClusterCollection> &lCCH,
        const edm::Handle<SimClusterCollection> &sCCH,
        const edm::Handle<CaloParticleCollection> &cPCH,
        const edm::Handle<std::map<uint, std::vector<uint>>> &simTrackstersMapH,
        const edm::Handle<ticl::TracksterCollection> &sTCH,
        const edm::Handle<ticl::TracksterCollection> &sTfromCPCH,
        const validationType valType) const;
  };
}  // namespace hgcal

#endif
