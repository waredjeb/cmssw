#ifndef SimDataFormats_Associations_TracksterToSimTracksterHitLCAssociator_h
#define SimDataFormats_Associations_TracksterToSimTracksterHitLCAssociator_h
// Original Author:  Leonardo Cristella

// system include files
#include <memory>

// user include files

#include "SimDataFormats/Associations/interface/TracksterToSimTracksterHitLCAssociatorBaseImpl.h"

// forward declarations

namespace hgcal {

  class TracksterToSimTracksterHitLCAssociator {
  public:
    TracksterToSimTracksterHitLCAssociator(std::unique_ptr<hgcal::TracksterToSimTracksterHitLCAssociatorBaseImpl>);
    TracksterToSimTracksterHitLCAssociator() = default;
    TracksterToSimTracksterHitLCAssociator(TracksterToSimTracksterHitLCAssociator &&) = default;
    TracksterToSimTracksterHitLCAssociator &operator=(TracksterToSimTracksterHitLCAssociator &&) = default;
    TracksterToSimTracksterHitLCAssociator(const TracksterToSimTracksterHitLCAssociator &) = delete;  // stop default
    const TracksterToSimTracksterHitLCAssociator &operator=(const TracksterToSimTracksterHitLCAssociator &) =
        delete;  // stop default

    ~TracksterToSimTracksterHitLCAssociator() = default;

    // ---------- const member functions ---------------------
    /// Associate a Trackster to SimClusters
    hgcal::RecoToSimCollectionSimTracksters associateRecoToSim(
        const edm::Handle<ticl::TracksterCollection> &tCH,
        const edm::Handle<reco::CaloClusterCollection> &lCCH,
        const edm::Handle<SimClusterCollection> &sCCH,
        const edm::Handle<CaloParticleCollection> &cPCH,
        const edm::Handle<std::map<uint, std::vector<uint>>> &simTrackstersMapH,
        const edm::Handle<ticl::TracksterCollection> &sTCH,
        const edm::Handle<ticl::TracksterCollection> &sTfromCPCH,
        const validationType valType) const {
      return m_impl->associateRecoToSim(tCH, lCCH, sCCH, cPCH, simTrackstersMapH, sTCH, sTfromCPCH, valType);
    };

    /// Associate a SimCluster to Tracksters
    hgcal::SimToRecoCollectionSimTracksters associateSimToReco(
        const edm::Handle<ticl::TracksterCollection> &tCH,
        const edm::Handle<reco::CaloClusterCollection> &lCCH,
        const edm::Handle<SimClusterCollection> &sCCH,
        const edm::Handle<CaloParticleCollection> &cPCH,
        const edm::Handle<std::map<uint, std::vector<uint>>> &simTrackstersMapH,
        const edm::Handle<ticl::TracksterCollection> &sTCH,
        const edm::Handle<ticl::TracksterCollection> &sTfromCPCH,
        const validationType valType) const {
      return m_impl->associateSimToReco(tCH, lCCH, sCCH, cPCH, simTrackstersMapH, sTCH, sTfromCPCH, valType);
    }

  private:
    // ---------- member data --------------------------------
    std::unique_ptr<TracksterToSimTracksterHitLCAssociatorBaseImpl> m_impl;
  };
}  // namespace hgcal

#endif
