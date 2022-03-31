// Original Author: Leonardo Cristella

#include "SimDataFormats/Associations/interface/TracksterToSimTracksterHitLCAssociatorBaseImpl.h"

namespace hgcal {
  TracksterToSimTracksterHitLCAssociatorBaseImpl::TracksterToSimTracksterHitLCAssociatorBaseImpl(){};
  TracksterToSimTracksterHitLCAssociatorBaseImpl::~TracksterToSimTracksterHitLCAssociatorBaseImpl(){};

  hgcal::RecoToSimCollectionSimTracksters TracksterToSimTracksterHitLCAssociatorBaseImpl::associateRecoToSim(
      const edm::Handle<ticl::TracksterCollection> &tCH,
      const edm::Handle<reco::CaloClusterCollection> &lCCH,
      const edm::Handle<SimClusterCollection> &sCCH,
      const edm::Handle<CaloParticleCollection> &cPCH,
      const edm::Handle<std::map<uint, std::vector<uint>>> &simTrackstersMapH,
      const edm::Handle<ticl::TracksterCollection> &sTCH,
      const edm::Handle<ticl::TracksterCollection> &sTfromCPCH,
      const validationType valType) const {
    return hgcal::RecoToSimCollectionSimTracksters();
  }

  hgcal::SimToRecoCollectionSimTracksters TracksterToSimTracksterHitLCAssociatorBaseImpl::associateSimToReco(
      const edm::Handle<ticl::TracksterCollection> &tCH,
      const edm::Handle<reco::CaloClusterCollection> &lCCH,
      const edm::Handle<SimClusterCollection> &sCCH,
      const edm::Handle<CaloParticleCollection> &cPCH,
      const edm::Handle<std::map<uint, std::vector<uint>>> &simTrackstersMapH,
      const edm::Handle<ticl::TracksterCollection> &sTCH,
      const edm::Handle<ticl::TracksterCollection> &sTfromCPCH,
      const validationType valType) const {
    return hgcal::SimToRecoCollectionSimTracksters();
  }

}  // namespace hgcal
