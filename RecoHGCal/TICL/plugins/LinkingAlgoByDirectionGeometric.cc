#include <cmath>
#include <string>
#include "RecoHGCal/TICL/plugins/LinkingAlgoByDirectionGeometric.h"

#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Common.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "TrackstersPCA.h"
using namespace ticl;

LinkingAlgoByDirectionGeometric::LinkingAlgoByDirectionGeometric(const edm::ParameterSet &conf)
    : LinkingAlgoBase(conf),
      del_tk_ts_layer1_(conf.getParameter<double>("delta_tk_ts_layer1")),
      del_tk_ts_int_(conf.getParameter<double>("delta_tk_ts_interface")),
      del_ts_em_had_(conf.getParameter<double>("delta_ts_em_had")),
      del_ts_had_had_(conf.getParameter<double>("delta_ts_had_had")),
      separationSmall_threshold_(conf.getParameter<double>("separationSmall")),
      separation_threshold_(conf.getParameter<double>("separation")),
      maxDepth_(conf.getParameter<int>("maxDepth")),
      timing_quality_threshold_(conf.getParameter<double>("track_time_quality_threshold")),
      cutTk_(conf.getParameter<std::string>("cutTk")) {}

LinkingAlgoByDirectionGeometric::~LinkingAlgoByDirectionGeometric() {}

void LinkingAlgoByDirectionGeometric::initialize(const HGCalDDDConstants *hgcons,
                                                 const hgcal::RecHitTools rhtools,
                                                 const edm::ESHandle<MagneticField> bfieldH,
                                                 const edm::ESHandle<Propagator> propH) {
  hgcons_ = hgcons;
  rhtools_ = rhtools;
  buildLayers();

  bfield_ = bfieldH;
  propagator_ = propH;
}

void LinkingAlgoByDirectionGeometric::fillTrackstersTile(const std::vector<Trackster> &tracksters,
                                                         std::array<TICLLayerTile, 2> &tracksterTiles) {
  // needs only the positive Z co-ordinate of the surface to propagate to
  // the correct sign is calculated inside according to the barycenter of trackster

  unsigned int idx = 0;
  for (const auto &t : tracksters) {
    Vector const &baryc = t.barycenter();

    if (baryc.Eta() > 0)
      tracksterTiles[1].fill(baryc.eta(), baryc.phi(), idx);

    else if (baryc.Eta() < 0)
      tracksterTiles[0].fill(baryc.eta(), baryc.phi(), idx);
    idx++;
  }
}

void fillTrackstersSmallTile(const std::vector<Trackster> &tracksters,
                             const std::vector<int> &maskSmall,
                             std::array<TICLLayerTile, 2> &tracksterTiles) {
  // needs only the positive Z co-ordinate of the surface to propagate to
  // the correct sign is calculated inside according to the barycenter of trackster

  for (size_t idx = 0; idx < tracksters.size(); idx++) {
    if (maskSmall[idx]) {
      auto const &t = tracksters[idx];
      auto const &baryc = t.barycenter();

      if (baryc.Eta() > 0)
        tracksterTiles[1].fill(baryc.eta(), baryc.phi(), idx);

      else if (baryc.Eta() < 0)
        tracksterTiles[0].fill(baryc.eta(), baryc.phi(), idx);
    }
  }
}
void LinkingAlgoByDirectionGeometric::fillTracksTile(const std::vector<std::pair<Vector, unsigned>> &propTracks,
                                                     std::array<TICLLayerTile, 2> &tracksTiles) {
  // needs only the positive Z co-ordinate of the surface to propagate to
  // the correct sign is calculated inside according to the barycenter of trackster

  for (const auto &t : propTracks) {
    Vector const &vec = t.first;

    if (vec.Eta() > 0)
      tracksTiles[1].fill(vec.eta(), vec.phi(), t.second);

    else if (vec.Eta() < 0)
      tracksTiles[0].fill(vec.eta(), vec.phi(), t.second);
  }
}
void LinkingAlgoByDirectionGeometric::findSmallTrackstersInWindow(const std::vector<Trackster> &tracksters,
                                                                  const std::vector<int> &maskSmall,
                                                                  const std::array<TICLLayerTile, 2> &tracksterTiles,
                                                                  const float delta,
                                                                  const float separation,
                                                                  std::vector<std::vector<unsigned>> &resultCollection,
                                                                  std::vector<float> &distancesVec,
                                                                  std::vector<int> &distancesVecIdx,
                                                                  bool useMask = false) {
  std::vector<int> mask(tracksters.size(), 0);
  float delta2 = delta * delta;
  float separation2 = separation * separation;
  auto distance = [](float r0, float r1, float phi0, float phi1) {
    auto delta_phi = reco::deltaPhi(phi0, phi1);
    return (r0 - r1) * (r0 - r1) + r1 * r1 * delta_phi * delta_phi;
  };
  for (size_t tracksterIndex = 0; tracksterIndex < tracksters.size(); tracksterIndex++) {
    if (!maskSmall[tracksterIndex]) {  //big tracksters
      auto const &t = tracksters[tracksterIndex];
      auto const &barycenter = t.barycenter();
      float trackster_eta = barycenter.eta();
      float trackster_phi = barycenter.phi();
      float trackster_z = barycenter.z();
      float trackster_r_over_absz =
          sqrt(barycenter.x() * barycenter.x() + barycenter.y() * barycenter.y()) / std::abs(trackster_z);
      auto sideZ = trackster_eta > 0;  //forward or backward region
      const TICLLayerTile &tile = tracksterTiles[sideZ];
      float eta_min = std::max(abs(trackster_eta) - delta, (float)TileConstants::minEta);
      float eta_max = std::min(abs(trackster_eta) + delta, (float)TileConstants::maxEta);
      // get range of bins touched by delta
      std::array<int, 4> search_box =
          tile.searchBoxEtaPhi(eta_min, eta_max, trackster_phi - delta, trackster_phi + delta);
      std::vector<unsigned> in_delta;
      std::vector<float> distances2;
      int trackster_n = 1;
      for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
        for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
          const auto &in_tile = tile[tile.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];
          for (const unsigned &t_i : in_tile) {
            if (t_i != tracksterIndex) {
              // calculate actual distances of tracksters to the seed for a more accurate cut
              auto const &tracksterInTile = tracksters[t_i];
              auto const &tracksterInTile_barycenter = tracksterInTile.barycenter();
              float tracksterInTile_eta = tracksterInTile_barycenter.eta();
              float tracksterInTile_phi = tracksterInTile_barycenter.phi();
              float tracksterInTile_z = tracksterInTile_barycenter.z();
              float tracksterInTile_r_over_absz =
                  sqrt(tracksterInTile_barycenter.x() * tracksterInTile_barycenter.x() +
                       tracksterInTile_barycenter.y() * tracksterInTile_barycenter.y()) /
                  std::abs(tracksterInTile_z);
              //            //std::cout << "Trackster Name\tTrackster X\tTrackster Y\tTrackster Z\tTrackster ROverAbsZ\t Trackster R "
              //                      << std::endl;
              //            //std::cout << "Trackster Ref \t" << barycenter.x() << "\t" << barycenter.y() << "\t" << barycenter.z()
              //                      << "\t" << trackster_r_over_absz << "\t" << trackster_r_over_absz * trackster_z << std::endl;
              //
              //            //std::cout << "Trackster " << trackster_n << "\t" << tracksterInTile_barycenter.x() << "\t"
              //                      << tracksterInTile_barycenter.y() << "\t" << tracksterInTile_barycenter.z() << "\t"
              //                      << tracksterInTile_r_over_absz << "\t" << tracksterInTile_r_over_absz * tracksterInTile_z
              //                      << std::endl;
              auto sep2 = distance(trackster_r_over_absz * trackster_z,
                                   tracksterInTile_r_over_absz * trackster_z,
                                   trackster_phi,
                                   tracksterInTile_phi);

              trackster_n++;
              distancesVec.push_back(sep2);
              distancesVecIdx.push_back(tracksterIndex);
              if (sep2 < separation2) {
                in_delta.push_back(t_i);
                distances2.push_back(sep2);
              }
            }
          }
        }
      }
      // sort tracksters found in ascending order of their distances from the seed
      std::vector<unsigned> indices(in_delta.size());
      std::iota(indices.begin(), indices.end(), 0);
      std::sort(indices.begin(), indices.end(), [&](unsigned i, unsigned j) { return distances2[i] < distances2[j]; });

      // push back sorted tracksters in the result collection
      for (const unsigned &index : indices) {
        const auto &t_i = in_delta[index];
        if (!mask[t_i]) {
          resultCollection[tracksterIndex].push_back(t_i);
          if (useMask)
            mask[t_i] = 1;
        }
      }
    }
  }
  //std::cout << "Result Collection size " << resultCollection.size() << std::endl;
}
void LinkingAlgoByDirectionGeometric::findTrackstersInWindow(const std::vector<Trackster> &tracksters,
                                                             const std::array<TICLLayerTile, 2> &tracksterTiles,
                                                             const float delta,
                                                             const float separation,
                                                             std::vector<std::vector<unsigned>> &resultCollection,
                                                             std::vector<float> &distancesVec,
                                                             std::vector<int> &distancesVecIdx,
                                                             bool useMask = false) {
  std::vector<int> mask(tracksters.size(), 0);
  float delta2 = delta * delta;
  float separation2 = separation * separation;
  auto distance = [](float r0, float r1, float phi0, float phi1) {
    auto delta_phi = reco::deltaPhi(phi0, phi1);
    return (r0 - r1) * (r0 - r1) + r1 * r1 * delta_phi * delta_phi;
  };
  unsigned int seedId = 0;
  int tracksterIndex = 0;
  for (auto const &t : tracksters) {
    auto const &barycenter = t.barycenter();
    float trackster_eta = barycenter.eta();
    float trackster_phi = barycenter.phi();
    float trackster_z = barycenter.z();
    float trackster_r_over_absz =
        sqrt(barycenter.x() * barycenter.x() + barycenter.y() * barycenter.y()) / std::abs(trackster_z);
    auto sideZ = trackster_eta > 0;  //forward or backward region
    const TICLLayerTile &tile = tracksterTiles[sideZ];
    float eta_min = std::max(abs(trackster_eta) - delta, (float)TileConstants::minEta);
    float eta_max = std::min(abs(trackster_eta) + delta, (float)TileConstants::maxEta);
    // get range of bins touched by delta
    std::array<int, 4> search_box =
        tile.searchBoxEtaPhi(eta_min, eta_max, trackster_phi - delta, trackster_phi + delta);
    std::vector<unsigned> in_delta;
    std::vector<float> distances2;
    int trackster_n = 1;
    for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
      for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
        const auto &in_tile = tile[tile.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];
        for (const unsigned &t_i : in_tile) {
          if (t_i != seedId) {
            // calculate actual distances of tracksters to the seed for a more accurate cut
            auto const &tracksterInTile = tracksters[t_i];
            auto const &tracksterInTile_barycenter = tracksterInTile.barycenter();
            float tracksterInTile_eta = tracksterInTile_barycenter.eta();
            float tracksterInTile_phi = tracksterInTile_barycenter.phi();
            float tracksterInTile_z = tracksterInTile_barycenter.z();
            float tracksterInTile_r_over_absz = sqrt(tracksterInTile_barycenter.x() * tracksterInTile_barycenter.x() +
                                                     tracksterInTile_barycenter.y() * tracksterInTile_barycenter.y()) /
                                                std::abs(tracksterInTile_z);
            //            //std::cout << "Trackster Name\tTrackster X\tTrackster Y\tTrackster Z\tTrackster ROverAbsZ\t Trackster R "
            //                      << std::endl;
            //            //std::cout << "Trackster Ref \t" << barycenter.x() << "\t" << barycenter.y() << "\t" << barycenter.z()
            //                      << "\t" << trackster_r_over_absz << "\t" << trackster_r_over_absz * trackster_z << std::endl;
            //
            //            //std::cout << "Trackster " << trackster_n << "\t" << tracksterInTile_barycenter.x() << "\t"
            //                      << tracksterInTile_barycenter.y() << "\t" << tracksterInTile_barycenter.z() << "\t"
            //                      << tracksterInTile_r_over_absz << "\t" << tracksterInTile_r_over_absz * tracksterInTile_z
            //                      << std::endl;
            auto sep2 = distance(trackster_r_over_absz * trackster_z,
                                 tracksterInTile_r_over_absz * trackster_z,
                                 trackster_phi,
                                 tracksterInTile_phi);

            trackster_n++;
            distancesVec.push_back(sep2);
            distancesVecIdx.push_back(tracksterIndex);
            if (sep2 < separation2) {
              in_delta.push_back(t_i);
              distances2.push_back(sep2);
            }
          }
        }
      }
    }
    // sort tracksters found in ascending order of their distances from the seed
    std::vector<unsigned> indices(in_delta.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](unsigned i, unsigned j) { return distances2[i] < distances2[j]; });

    // push back sorted tracksters in the result collection
    for (const unsigned &index : indices) {
      const auto &t_i = in_delta[index];
      if (!mask[t_i]) {
        resultCollection[seedId].push_back(t_i);
        if (useMask)
          mask[t_i] = 1;
      }
    }
    seedId++;
    tracksterIndex++;
  }
  //std::cout << "Result Collection size " << resultCollection.size() << std::endl;
}

void LinkingAlgoByDirectionGeometric::tracksterToTrackLinking(
    std::vector<Trackster> &tracksterMergeCollection,
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
    std::vector<TICLCandidate> &candidates,
    std::vector<TICLCandidate> &chargedCandidatesFromTracks,
    std::vector<float> &separations2_for_ntuples,
    std::vector<float> &separations2_ETCompatible_for_ntuples) {
  std::vector<int> mask(tracksterMergeTiles.size(), 1);
  std::vector<int> maskMergeCollection(tracksterMergeCollection.size(), 1);
  std::vector<int> maskTracks(propTracks.size(), 1);
  std::vector<Trackster> orderedTracksterMergeCollection;
  float delta2 = delta * delta;
  float separation2 = separation * separation;
  auto distance = [](float r0, float r1, float phi0, float phi1) {
    auto delta_phi = reco::deltaPhi(phi0, phi1);
    return (r0 - r1) * (r0 - r1) + r1 * r1 * delta_phi * delta_phi;
  };
  unsigned int seedId = 0;

  for (size_t track_i; track_i < propTracks.size(); track_i++) {
    auto const &propTrack = propTracks[track_i];

    float bestSep2 = 9999;
    int bestMergeTracksterIndex = -1;  //initialize as invalid index
    auto const &trackPosition = propTrack.first;

    float track_eta = trackPosition.eta();

    float track_phi = trackPosition.phi();

    float track_z = trackPosition.z();

    float track_r_over_absz =
        sqrt(trackPosition.x() * trackPosition.x() + trackPosition.y() * trackPosition.y()) / std::abs(track_z);
    auto sideZ = track_eta > 0;  //forward or backward region

    const TICLLayerTile &tile = tracksterMergeTiles[sideZ];

    float eta_min = std::max(abs(track_eta) - delta, (float)TileConstants::minEta);
    float eta_max = std::min(abs(track_eta) + delta, (float)TileConstants::maxEta);
    // get range of bins touched by delta

    std::array<int, 4> search_box = tile.searchBoxEtaPhi(eta_min, eta_max, track_phi - delta, track_phi + delta);

    std::vector<unsigned> in_delta;
    std::vector<float> distances2;
    int track_n = 1;
    for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
      for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
        const auto &in_tile = tile[tile.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];
        for (const unsigned &t_i : in_tile) {
          // calculate actual distances of tracksters to the seed for a more accurate cut
          auto const &tracksterInTile = tracksterMergeCollection[t_i];
          auto const &tracksterInTile_barycenter = tracksterInTile.barycenter();
          float tracksterInTile_eta = tracksterInTile_barycenter.eta();
          float tracksterInTile_phi = tracksterInTile_barycenter.phi();
          float tracksterInTile_z = tracksterInTile_barycenter.z();
          float tracksterInTile_r_over_absz = sqrt(tracksterInTile_barycenter.x() * tracksterInTile_barycenter.x() +
                                                   tracksterInTile_barycenter.y() * tracksterInTile_barycenter.y()) /
                                              std::abs(tracksterInTile_z);
          //          //std::cout << "Trackster Name\tTrackster X\tTrackster Y\tTrackster Z\tTrackster ROverAbsZ\t Trackster R "
          //                    << std::endl;
          //          //std::cout << "Track Ref \t" << trackPosition.x() << "\t" << trackPosition.y() << "\t" << trackPosition.z()
          //                    << "\t" << track_r_over_absz << "\t" << track_r_over_absz * tracksterInTile_z << std::endl;
          //
          //          //std::cout << "Trackster " << track_n << "\t" << tracksterInTile_barycenter.x() << "\t"
          //                    << tracksterInTile_barycenter.y() << "\t" << tracksterInTile_barycenter.z() << "\t"
          //                    << tracksterInTile_r_over_absz << "\t" << tracksterInTile_r_over_absz * tracksterInTile_z
          //                    << std::endl;
          auto sep2 = distance(
              track_r_over_absz * track_z, tracksterInTile_r_over_absz * track_z, track_phi, tracksterInTile_phi);

          track_n++;
          auto tkRef = reco::TrackRef(tkH, propTrack.second);
          auto track_time = tkTime[tkRef];
          auto track_timeErr = tkTimeErr[tkRef];
          auto track_timeQual = tkTimeQual[tkRef];
          auto energyTimeCompatible = timeAndEnergyCompatible(
              tracks[propTrack.second], tracksterMergeCollection[t_i], track_time, track_timeErr, track_timeQual);
          if (!energyTimeCompatible) {
            separations2_for_ntuples.push_back(sep2);
          } else {
            separations2_ETCompatible_for_ntuples.push_back(sep2);
          }
          if (sep2 < separation2 && sep2 < bestSep2 && energyTimeCompatible) {
            bestSep2 = sep2;
            bestMergeTracksterIndex = t_i;
            distances2.push_back(sep2);
          }
        }
      }
    }

    // sort tracksters found in ascending order of their distances from the seed
    if (bestMergeTracksterIndex != -1) {  //matched track with trackster --> build charged candidate

      if (maskTracks[track_i] && maskMergeCollection[bestMergeTracksterIndex]) {
        TICLCandidate chargedCandidate;

        for (auto const &t_indices : tracksterMergeCollectionIndices[bestMergeTracksterIndex]) {
          chargedCandidate.addTrackster(edm::Ptr<Trackster>(tsH, t_indices));
        }
        chargedCandidate.setTrackPtr(edm::Ptr<reco::Track>(tkH, propTrack.second));
        candidates.push_back(chargedCandidate);
        orderedTracksterMergeCollection.push_back(tracksterMergeCollection[bestMergeTracksterIndex]);
      }

      maskMergeCollection[bestMergeTracksterIndex] = 0;

      maskTracks[track_i] = 0;
      auto bestTrackster = tracksterMergeCollection[bestMergeTracksterIndex];

      //std::cout << "Track Ref Energy " << tracks[track_i].p() << " Position " << propTrack.first.x() << " "
      //       << propTrack.first.y() << " " << propTrack.first.z() << std::endl;
      //std::cout << "Matched Trackster Energy " << bestTrackster.raw_energy() << " Position "
      //      << bestTrackster.barycenter().x() << " " << bestTrackster.barycenter().y() << " "
      //   << bestTrackster.barycenter().z() << std::endl;
    }
    seedId++;
  }
  //promote all the non used tracks to charged candidates
  for (size_t mask_track_i = 0; mask_track_i < maskTracks.size(); mask_track_i++) {
    if (maskTracks[mask_track_i]) {
      TICLCandidate chargedCandidateFromTk;
      chargedCandidateFromTk.setTrackPtr(edm::Ptr<reco::Track>(tkH, propTracks[mask_track_i].second));
      chargedCandidatesFromTracks.push_back(chargedCandidateFromTk);
    }
  }

  for (size_t mask_trackster_i = 0; mask_trackster_i < maskMergeCollection.size(); mask_trackster_i++) {
    if (maskMergeCollection[mask_trackster_i]) {
      TICLCandidate neutralCandidate;
      auto const &tracksterIndices = tracksterMergeCollectionIndices[mask_trackster_i];
      for (auto const &tracksterIndex : tracksterIndices) {
        neutralCandidate.addTrackster(edm::Ptr<Trackster>(tsH, tracksterIndex));
      }
      candidates.push_back(neutralCandidate);
      orderedTracksterMergeCollection.push_back(tracksterMergeCollection[mask_trackster_i]);
    }
  }
  tracksterMergeCollection = orderedTracksterMergeCollection;
}

bool LinkingAlgoByDirectionGeometric::timeCompatible(const Trackster &trackster1, const Trackster &trackster2) {
  float tsT1 = trackster1.time();
  float tsTErr1 = trackster1.timeError();

  float tsT2 = trackster2.time();
  float tsTErr2 = trackster2.timeError();
  auto timeCompatible = false;
  if (tsT1 != -99 && tsT2 != -99) {
    timeCompatible = (std::abs(tsT1 - tsT2) < maxDeltaT_ * sqrt(tsTErr1 * tsTErr1 + tsTErr2 * tsTErr2));
  } else {
    timeCompatible = true;
  }
  return timeCompatible;
}
bool LinkingAlgoByDirectionGeometric::timeAndEnergyCompatible(const reco::Track &track,
                                                              const Trackster &trackster,
                                                              const float &tkT,
                                                              const float &tkTErr,
                                                              const float &tkTimeQual) {
  float threshold = std::min(0.2 * trackster.regressed_energy(), 10.0);

  bool energyCompatible = (trackster.regressed_energy() < track.p() + threshold);
  // compatible if trackster time is within 3sigma of
  // track time; compatible if either: no time assigned
  // to trackster or track time quality is below threshold
  float tsT = trackster.time();
  float tsTErr = trackster.timeError();

  bool timeCompatible = false;

  if (tsT == -99. or tkTimeQual < timing_quality_threshold_)
    timeCompatible = true;
  else {
    timeCompatible = (std::abs(tsT - tkT) < maxDeltaT_ * sqrt(tsTErr * tsTErr + tkTErr * tkTErr));
  }

  if (LinkingAlgoBase::algo_verbosity_ > VerbosityLevel::Advanced) {
    if (!(energyCompatible))
      LogDebug("LinkingAlgoByDirectionGeometric") << "energy incompatible : track p " << track.p()
                                                  << " trackster energy " << trackster.regressed_energy() << "\n";
    if (!(timeCompatible))
      LogDebug("LinkingAlgoByDirectionGeometric") << "time incompatible : track time " << tkT << " +/- " << tkTErr
                                                  << " trackster time " << tsT << " +/- " << tsTErr << "\n";
  }
  return energyCompatible && timeCompatible;
}

void LinkingAlgoByDirectionGeometric::recordTrackster(const unsigned ts,  //trackster index
                                                      const std::vector<Trackster> &tracksters,
                                                      const edm::Handle<std::vector<Trackster>> tsH,
                                                      std::vector<unsigned> &ts_mask,
                                                      float &energy_in_candidate,
                                                      TICLCandidate &candidate) {}

void LinkingAlgoByDirectionGeometric::buildLayers() {
  // build disks at HGCal front & EM-Had interface for track propagation

  float zVal = hgcons_->waferZ(1, true);
  std::pair<float, float> rMinMax = hgcons_->rangeR(zVal, true);

  float zVal_interface = rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z();
  std::pair<float, float> rMinMax_interface = hgcons_->rangeR(zVal_interface, true);

  for (int iSide = 0; iSide < 2; ++iSide) {
    float zSide = (iSide == 0) ? (-1. * zVal) : zVal;
    firstDisk_[iSide] =
        std::make_unique<GeomDet>(Disk::build(Disk::PositionType(0, 0, zSide),
                                              Disk::RotationType(),
                                              SimpleDiskBounds(rMinMax.first, rMinMax.second, zSide - 0.5, zSide + 0.5))
                                      .get());

    zSide = (iSide == 0) ? (-1. * zVal_interface) : zVal_interface;
    interfaceDisk_[iSide] = std::make_unique<GeomDet>(
        Disk::build(Disk::PositionType(0, 0, zSide),
                    Disk::RotationType(),
                    SimpleDiskBounds(rMinMax_interface.first, rMinMax_interface.second, zSide - 0.5, zSide + 0.5))
            .get());
  }
}
void LinkingAlgoByDirectionGeometric::dumpLinksFound(std::vector<std::vector<unsigned>> &resultCollection,
                                                     const char *label) const {
  //#ifdef EDM_ML_DEBUG
  //  if (!(LinkingAlgoBase::algo_verbosity_ > VerbosityLevel::Advanced))
  //    return;

  //std::cout << "All links found - " << label << "\n";
  //std::cout << "(seed can either be a track or trackster depending on the step)\n";
  for (unsigned i = 0; i < resultCollection.size(); ++i) {
    //std::cout << "seed " << i << " - tracksters : ";
    const auto &links = resultCollection[i];
    for (unsigned j = 0; j < links.size(); ++j) {
      //std::cout << j;
    }
    //std::cout << "\n";
  }
  //#endif  // EDM_ML_DEBUG
}

void DFSVisits(std::vector<std::vector<unsigned>> &graph,
               std::vector<bool> &visits,
               Trackster &outTrackster,
               std::vector<unsigned> &outTracksterIndices,
               std::vector<std::vector<unsigned>> &resultCollectionIndices,
               const std::vector<Trackster> &tracksters,
               const std::vector<reco::CaloCluster> &layerClusters,
               int u,
               std::string &tabs,
               int &currentDepth,
               const int maxDepth,
               bool secondStep) {
  //  //std::cout << "Visiting " << u << std::endl;
  auto updatedSize = outTrackster.vertices().size();
  //  outTracksterIndices.push_back(u);
  for (auto j = graph[u].begin(); j != graph[u].end(); j++) {
    if (!visits.at(*j) && currentDepth <= maxDepth) {
      currentDepth += 1;
      visits[*j] = true;
      auto const &thisTrackster = tracksters[*j];

      tabs.push_back('\t');
      updatedSize += thisTrackster.vertices().size();
      outTrackster.vertices().reserve(updatedSize);
      outTrackster.vertex_multiplicity().reserve(updatedSize);
      std::copy(std::begin(thisTrackster.vertices()),
                std::end(thisTrackster.vertices()),
                std::back_inserter(outTrackster.vertices()));
      std::copy(std::begin(thisTrackster.vertex_multiplicity()),
                std::end(thisTrackster.vertex_multiplicity()),
                std::back_inserter(outTrackster.vertex_multiplicity()));
      //std::cout << tabs << "Visiting " << *j << " Energy " << thisTrackster.raw_energy() << " Position "
        //        << thisTrackster.barycenter() << std::endl;

      if (!secondStep) {
        outTracksterIndices.push_back(*j);
      } else {
        for (auto const &v : resultCollectionIndices[*j]) {
          outTracksterIndices.push_back(v);
        }
      }
      DFSVisits(graph,
                visits,
                outTrackster,
                outTracksterIndices,
                resultCollectionIndices,
                tracksters,
                layerClusters,
                *j,
                tabs,
                currentDepth,
                maxDepth,
                secondStep);
      if (j == graph[u].end() - 1) {
        tabs.pop_back();
      }
    }
  }
}
// DFS to collect all the links
void DFS(std::vector<std::vector<unsigned>> &graph,
         const std::vector<Trackster> &tracksters,
         const std::vector<reco::CaloCluster> &layerClusters,
         std::vector<Trackster> &resultCollection,
         std::vector<std::vector<unsigned>> &resultCollectionIndices,
         const int maxDepth,
         bool secondStep) {
  int graphSize = graph.size();
  std::vector<bool> visits(graphSize, false);
  //std::cout << "Initial Graph size " << graphSize << " Visits size " << visits.size() << std::endl;
  for (int i = 0; i < graphSize; i++) {
    std::string tabs;
    if (!visits[i]) {
      visits[i] = true;
      //      resultCollectionIndices[i].push_back(i);
      tabs.push_back('\t');
      Trackster outTrackster = tracksters[i];
      std::vector<unsigned> outTracksterIndices = {static_cast<unsigned>(i)};

      int currentDepth = 0;
      //std::cout << " -- Trackster " << i << " Energy " << outTrackster.raw_energy() << " Position "
        //        << outTrackster.barycenter() << std::endl;
      DFSVisits(graph,
                visits,
                outTrackster,
                outTracksterIndices,
                resultCollectionIndices,
                tracksters,
                layerClusters,
                i,
                tabs,
                currentDepth,
                maxDepth,
                secondStep);

      resultCollection.push_back(outTrackster);
      resultCollectionIndices.push_back(outTracksterIndices);
    }
  }
  resultCollection.shrink_to_fit();
}

void DFSFinal(std::vector<std::vector<unsigned>> &graph,
              const std::vector<Trackster> &tracksters,
              const std::vector<reco::CaloCluster> &layerClusters,
              std::vector<std::vector<unsigned>> &previousCollectionIndices,
              std::vector<Trackster> &resultCollection,
              std::vector<std::vector<unsigned>> &resultCollectionIndices,
              const int maxDepth,
              bool secondStep) {
  int graphSize = graph.size();
  std::vector<bool> visits(graphSize, false);
  for (int i = 0; i < graphSize; i++) {
    std::string tabs;
    if (!visits[i]) {
      visits[i] = true;
      //      resultCollectionIndices[i].push_back(i);
      tabs.push_back('\t');
      Trackster outTrackster = tracksters[i];
      std::vector<unsigned> outTracksterIndices = previousCollectionIndices[i];
      //std::cout << " -- Trackster " << i << " Energy " << outTrackster.raw_energy() << " Position "
        //        << outTrackster.barycenter() << std::endl;
      auto currentDepth = 0;
      DFSVisits(graph,
                visits,
                outTrackster,
                outTracksterIndices,
                previousCollectionIndices,
                tracksters,
                layerClusters,
                i,
                tabs,
                currentDepth,
                maxDepth,
                secondStep);

      resultCollection.push_back(outTrackster);
      resultCollectionIndices.push_back(outTracksterIndices);
    }
  }
  resultCollection.shrink_to_fit();
  //std::cout << "RESULT COLLECTION SIZE " << resultCollection.size() << " RESULT COLLECTION INDICES SIZE " << resultCollectionIndices.size() << std::endl;
}

void findSmallTracksters(const std::vector<Trackster> &tracksters,
                         std::vector<int> &maskSmall,
                         const std::vector<reco::CaloCluster> &lcs,
                         const hgcal::RecHitTools &rhtools,
                         const int maxLayers) {
  for (size_t i = 0; i < tracksters.size(); i++) {
    std::vector<int> layers;
    auto const &t = tracksters[i];
    auto const &vertices = t.vertices();
    for (auto const &v : vertices) {
      auto const &hitId = lcs[v].hitsAndFractions()[0].first;
      auto const &layer = rhtools.getLayerWithOffset(hitId);
      layers.push_back(layer);
    }
    auto it = std::unique(layers.begin(), layers.end());
    layers.resize(std::distance(layers.begin(), it));
    maskSmall[i] = static_cast<int>(layers.size()) < maxLayers ? 1 : 0;
  }
}

void LinkingAlgoByDirectionGeometric::linkTracksters(const edm::Handle<std::vector<reco::Track>> tkH,
                                                     const edm::ValueMap<float> &tkTime,
                                                     const edm::ValueMap<float> &tkTimeErr,
                                                     const edm::ValueMap<float> &tkTimeQual,
                                                     const std::vector<reco::Muon> &muons,
                                                     const edm::Handle<std::vector<Trackster>> tsH,
                                                     const std::vector<reco::CaloCluster> &layerClusters,
                                                     const edm::ValueMap<std::pair<float, float>> &layerClustersTimes,
                                                     std::vector<Trackster> &resultTrackstersMerged,
                                                     std::vector<TICLCandidate> &candidates,
                                                     std::vector<TICLCandidate> &chargedCandidatesFromTracks,
                                                     const EnergyRegressionAndIDModel &model,
                                                     std::vector<float> &separations2_for_ntuples,
                                                     std::vector<float> &separations2_ETCompatible_for_ntuples,
                                                     std::vector<float> &distancesVec,
                                                     std::vector<int> &distancesVecIdx,
                                                     std::vector<double>& prop_tracks_x,
                                                     std::vector<double>& prop_tracks_y,
                                                     std::vector<double>& prop_tracks_z,
                                                     std::vector<double>& prop_tracks_eta,
                                                     std::vector<double>& prop_tracks_phi,
                                                     std::vector<double>& prop_tracks_px,
                                                     std::vector<double>& prop_tracks_py,
                                                     std::vector<double>& prop_tracks_pz,
                                                     std::vector<bool>& masked_tracks
                                                     ) {
  const auto &tracks = *tkH;
  const auto &tracksters = *tsH;

  auto bFieldProd = bfield_.product();
  const Propagator &prop = (*propagator_);
  //needed to keep track of trackster indices when building TICLCandidate
  std::vector<std::vector<unsigned>> tracksterMergeSmallCollectionIndices;
  std::vector<std::vector<unsigned>> tracksterMergeCollectionIndices;
  //needed to keep track of trackster indices when building TICLCandidate
  std::vector<Trackster> tracksterMergeSmallCollection;
  std::vector<Trackster> tracksterMergeCollection;
  //after cleaning needed to keep track of trackster indices when building TICLCandidate
  std::vector<Trackster> resultTrackstersSmallMerged;  //after cleaning.
  std::vector<std::vector<unsigned>> resultTrackstersSmallMergedIndices;
  std::vector<std::vector<unsigned>> resultTrackstersMergedIndices;

  std::array<TICLLayerTile, 2> tracksterPropTiles = {};       // all Tracksters
  std::array<TICLLayerTile, 2> tracksterSmallPropTiles = {};  // all Tracksters
  std::array<TICLLayerTile, 2> tracksPropTiles = {};          // all Tracks

  if (LinkingAlgoBase::algo_verbosity_ > VerbosityLevel::Advanced)
    LogDebug("LinkingAlgoByDirectionGeometric") << "------- Geometric Linking ------- \n";

  // Propagate tracks
  std::vector<unsigned> candidateTrackIds;
  candidateTrackIds.reserve(tracks.size());
  std::vector<std::pair<Vector, unsigned>> trackPColl;     // propagated track points and index of track in collection
  std::vector<std::pair<Vector, unsigned>> tkPropIntColl;  // tracks propagated to lastLayerEE
  trackPColl.reserve(tracks.size());
  tkPropIntColl.reserve(tracks.size());
  for (unsigned i = 0; i < tracks.size(); ++i) {
    const auto &tk = tracks[i];
    reco::TrackRef trackref = reco::TrackRef(tkH, i);

    // veto tracks associated to muons
    int muId = PFMuonAlgo::muAssocToTrack(trackref, muons);

    if (LinkingAlgoBase::algo_verbosity_ > VerbosityLevel::Advanced)
      LogDebug("LinkingAlgoByDirectionGeometric")
          << "track " << i << " - eta " << tk.eta() << " phi " << tk.phi() << " time " << tkTime[reco::TrackRef(tkH, i)]
          << " time qual " << tkTimeQual[reco::TrackRef(tkH, i)] << "  muid " << muId << "\n";

    if (!cutTk_((tk)) or muId != -1) {
      masked_tracks[i] = false;
      prop_tracks_x.push_back(999);
      prop_tracks_y.push_back(999);
      prop_tracks_z.push_back(999);
      prop_tracks_eta.push_back(999);
      prop_tracks_phi.push_back(999);
      prop_tracks_px.push_back(999);
      prop_tracks_py.push_back(999);
      prop_tracks_pz.push_back(999);
      continue;
    }
    // record tracks that can be used to make a ticlcandidate
    candidateTrackIds.push_back(i);

    // don't consider tracks below 2 GeV for linking
    if (std::sqrt(tk.p() * tk.p() + ticl::mpion2) < tkEnergyCut_) {
      masked_tracks[i] = false;
      prop_tracks_x.push_back(999);
      prop_tracks_y.push_back(999);
      prop_tracks_z.push_back(999);
      prop_tracks_eta.push_back(999);
      prop_tracks_phi.push_back(999);
      prop_tracks_px.push_back(999);
      prop_tracks_py.push_back(999);
      prop_tracks_pz.push_back(999);
      continue;
    }
    int iSide = int(tk.eta() > 0);
    const auto &fts = trajectoryStateTransform::outerFreeState((tk), bFieldProd);
    // to the HGCal front
    const auto &tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
    Vector trackPFirstLayer;
    float trackEnergy = sqrt(tracks[i].outerP() * tracks[i].outerP()  + mpion2);
    if (tsos.isValid()) {
      Vector trackP(tsos.globalPosition().x(), tsos.globalPosition().y(), tsos.globalPosition().z());
      trackPFirstLayer = trackP;
      prop_tracks_x.push_back(trackP.x());
      prop_tracks_y.push_back(trackP.y());
      prop_tracks_z.push_back(trackP.z());
      prop_tracks_eta.push_back(trackP.eta());
      prop_tracks_phi.push_back(trackP.phi());
      prop_tracks_px.push_back(tsos.globalMomentum().x());
      prop_tracks_py.push_back(tsos.globalMomentum().y());
      prop_tracks_pz.push_back(tsos.globalMomentum().z());
//      std::cout << "First Layer " << std::endl;
//      std::cout << "Track " << i << " Energy " << trackEnergy << " Global Momentum" << tsos.globalMomentum() << " Position " << tsos.globalPosition() << std::endl;
//      std::cout << "Track " << i << " Energy " << trackEnergy << " Global Direction " << tsos.globalDirection() << " Position " << tsos.globalPosition() << std::endl;
      masked_tracks[i] = true;
      trackPColl.emplace_back(trackP, i);
    } else {
      masked_tracks[i] = false;
      // std::cout << "PROPAGATION NOT VALID! " << std::endl;
    }
    // to lastLayerEE
    const auto &tsos_int = prop.propagate(fts, interfaceDisk_[iSide]->surface());
    Vector trackPIntLayer;
    if (tsos_int.isValid()) {
      Vector trackP(tsos_int.globalPosition().x(), tsos_int.globalPosition().y(), tsos_int.globalPosition().z());
      trackPIntLayer = trackP;
      tkPropIntColl.emplace_back(trackP, i);
     // prop_tracks_xInt.push_back(trackP.x());
     // prop_tracks_yInt.push_back(trackP.y());
     // prop_tracks_zInt.push_back(trackP.z());
     // prop_tracks_etaInt.push_back(trackP.eta());
     // prop_tracks_phiInt.push_back(trackP.phi());
     // prop_tracks_pxInt.push_back(tsos.globalMomentum().x());
     // prop_tracks_pyInt.push_back(tsos.globalMomentum().y());
     // prop_tracks_pzInt.push_back(tsos.globalMomentum().z());
     // masked_tracks[i] = true;
      trackPColl.emplace_back(trackP, i);
//      std::cout << "Interface Layer " << std::endl;
//      std::cout << "Track " << i << "Energy " << trackEnergy << " Global Momentum " << tsos_int.globalMomentum() << " Position " << tsos_int.globalPosition() << std::endl;
//      std::cout << "Track " << i << "Energy " << trackEnergy << " Global Direction " << tsos_int.globalDirection() << " Position " << tsos_int.globalPosition() << std::endl;
    } else {
      //masked_tracks[i] = false;
      // std::cout << "PROPAGATION NOT VALID! " << std::endl;
    }
//    Vector directionFromTwoPoints = trackPIntLayer - trackPFirstLayer;
 //   std::cout << "Direction from two points " << directionFromTwoPoints.Unit() << std::endl;
  }  // Tracks
  tkPropIntColl.shrink_to_fit();
  trackPColl.shrink_to_fit();
  candidateTrackIds.shrink_to_fit();
  // std::cout << "++++++++++++++++++++ SMALL MERGING +++++++++++++++++++++++++++++++++++++++++ " << std::endl;
  std::vector<int> maskSmall;
  maskSmall.resize(tracksters.size());
  findSmallTracksters(tracksters, maskSmall, layerClusters, rhtools_, 5);
  fillTrackstersSmallTile(tracksters, maskSmall, tracksterSmallPropTiles);
  std::vector<std::vector<unsigned>> trackstersResults(tracksters.size());
  findSmallTrackstersInWindow(tracksters,
                              maskSmall,
                              tracksterSmallPropTiles,
                              0.1,
                              separationSmall_threshold_,
                              trackstersResults,
                              distancesVec,
                              distancesVecIdx);
  dumpLinksFound(trackstersResults, "Tracksters Links");

  std::vector<int> maskTracksters(tracksters.size(), 1);

  //  tracksterMergeSmallCollectionIndices.resize(tracksters.size());
  DFS(trackstersResults,
      tracksters,
      layerClusters,
      tracksterMergeSmallCollection,
      tracksterMergeSmallCollectionIndices,
      maxDepth_,
      false);
  assert(tracksterMergeSmallCollectionIndices.size() == tracksterMergeSmallCollection.size());
  //  std::sort(tracksterMergeSmallCollectionIndices.begin(), tracksterMergeSmallCollectionIndices.end(),
  //            [&](const std::vector<unsigned>& a, const std::vector<unsigned> b) {
  //              return a.size() > b.size();
  //            });
  for (size_t k = 0; k < tracksterMergeSmallCollectionIndices.size(); k++) {
    //std::cout << "T " << k << " Indices ";
    for (auto const &v : tracksterMergeSmallCollectionIndices[k]) {
      // std::cout << v << " ";
    }
    //std::cout << "\n";
  }

  for (size_t i = 0; i < tracksterMergeSmallCollection.size(); i++) {
    auto &outTrackster = tracksterMergeSmallCollection[i];
    auto &orig_vtx = outTrackster.vertices();
    auto vtx_sorted{orig_vtx};
    std::sort(std::begin(vtx_sorted), std::end(vtx_sorted));
    for (unsigned int iLC = 1; iLC < vtx_sorted.size(); ++iLC) {
      if (vtx_sorted[iLC] == vtx_sorted[iLC - 1]) {
        // Clean up duplicate LCs
        const auto lcIdx = vtx_sorted[iLC];
        const auto firstEl = std::find(orig_vtx.begin(), orig_vtx.end(), lcIdx);
        const auto firstPos = std::distance(std::begin(orig_vtx), firstEl);
        auto iDup = std::find(std::next(firstEl), orig_vtx.end(), lcIdx);
        while (iDup != orig_vtx.end()) {
          orig_vtx.erase(iDup);
          outTrackster.vertex_multiplicity().erase(outTrackster.vertex_multiplicity().begin() +
                                                   std::distance(std::begin(orig_vtx), iDup));
          outTrackster.vertex_multiplicity()[firstPos] -= 1;
          iDup = std::find(std::next(firstEl), orig_vtx.end(), lcIdx);
        };
      }
    }

    outTrackster.zeroProbabilities();
    if (!outTrackster.vertices().empty()) {
      resultTrackstersSmallMerged.push_back(outTrackster);
    }
  }

  // Create TrackstersMerge with assignPCA
  // Create vector of Linked Tracksters
  // Pass this vector to assignPCA to have the final TrackstersMerged?
  resultTrackstersSmallMerged.shrink_to_fit();
  //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  //std::cout << "+++++++++++++++++++++++++ FINAL MERGING +++++++++++++++++++++++++++" << std::endl;
  fillTrackstersTile(resultTrackstersSmallMerged, tracksterPropTiles);
  std::vector<std::vector<unsigned>> trackstersFinalResults(resultTrackstersSmallMerged.size());
  findTrackstersInWindow(resultTrackstersSmallMerged,
                         tracksterPropTiles,
                         0.1,
                         separation_threshold_,
                         trackstersFinalResults,
                         distancesVec,
                         distancesVecIdx);
  //  std::cout << "Trackster size " << tracksters.size() << " Result Tracksters Small Merged "
  //     << resultTrackstersSmallMerged.size() << " Corresponding indices " << trackstersFinalResults.size()
  //     << std::endl;
  //  dumpLinksFound(trackstersFinalResults, "Tracksters Links");
  //  tracksterMergeCollectionIndices.resize(resultTrackstersSmallMerged.size());
  DFSFinal(trackstersFinalResults,
           resultTrackstersSmallMerged,
           layerClusters,
           tracksterMergeSmallCollectionIndices,
           tracksterMergeCollection,
           tracksterMergeCollectionIndices,
           maxDepth_,
           true);

  assert(tracksterMergeCollection.size() == tracksterMergeCollectionIndices.size());

  for (size_t k = 0; k < tracksterMergeCollectionIndices.size(); k++) {
    //std::cout << "T " << k << " Indices ";
    for (auto const &v : tracksterMergeCollectionIndices[k]) {
      //  std::cout << v << " ";
    }
    // std::cout << "\n";
  }
  //
  for (size_t i = 0; i < tracksterMergeCollection.size(); i++) {
    auto &outTrackster = tracksterMergeCollection[i];
    auto &outTracksterIndices = tracksterMergeCollectionIndices[i];
    auto &orig_vtx = outTrackster.vertices();
    auto vtx_sorted{orig_vtx};
    std::sort(std::begin(vtx_sorted), std::end(vtx_sorted));
    for (unsigned int iLC = 1; iLC < vtx_sorted.size(); ++iLC) {
      if (vtx_sorted[iLC] == vtx_sorted[iLC - 1]) {
        // Clean up duplicate LCs
        const auto lcIdx = vtx_sorted[iLC];
        const auto firstEl = std::find(orig_vtx.begin(), orig_vtx.end(), lcIdx);
        const auto firstPos = std::distance(std::begin(orig_vtx), firstEl);
        auto iDup = std::find(std::next(firstEl), orig_vtx.end(), lcIdx);
        while (iDup != orig_vtx.end()) {
          orig_vtx.erase(iDup);
          outTrackster.vertex_multiplicity().erase(outTrackster.vertex_multiplicity().begin() +
                                                   std::distance(std::begin(orig_vtx), iDup));
          outTrackster.vertex_multiplicity()[firstPos] -= 1;
          iDup = std::find(std::next(firstEl), orig_vtx.end(), lcIdx);
        };
      }
    }

    outTrackster.zeroProbabilities();
    if (!outTrackster.vertices().empty()) {
      resultTrackstersMerged.push_back(outTrackster);
      resultTrackstersMergedIndices.push_back(outTracksterIndices);
    }
  }
  assignPCAtoTracksters(
      resultTrackstersMerged, layerClusters, layerClustersTimes, rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z());
  model.energyRegressionAndID(layerClusters, resultTrackstersMerged);
  //std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  //Link with tracks and final candidate building
  //Refill Tiles with trackstersMerge
  std::array<TICLLayerTile, 2> tracksterMergePropTiles = {};  // all TrackstersMerge
  fillTrackstersTile(resultTrackstersMerged, tracksterMergePropTiles);
  fillTracksTile(trackPColl, tracksPropTiles);
  //  std::vector<TICLCandidate> candidates;
  //  std::vector<TICLCandidate> chargedCandidatesFromTracks;
  const float delta_tmp = 0.03;
  tracksterToTrackLinking(resultTrackstersMerged,
                          resultTrackstersMergedIndices,
                          tracksters,
                          tsH,
                          tracks,
                          tkH,
                          tkTime,
                          tkTimeErr,
                          tkTimeQual,
                          trackPColl,
                          tracksterMergePropTiles,
                          del_tk_ts_int_,
                          separation_threshold_,
                          candidates,
                          chargedCandidatesFromTracks,
                          separations2_for_ntuples,
                          separations2_ETCompatible_for_ntuples);
  //std::cout << "All Candidates " << candidates.size() << std::endl;
  //std::cout << "Charged Candidates From Tracks " << chargedCandidatesFromTracks.size() << std::endl;
}  // linkTracksters
void LinkingAlgoByDirectionGeometric::fillPSetDescription(edm::ParameterSetDescription &desc) {
  desc.add<std::string>("cutTk",
                        "1.48 < abs(eta) < 3.0 && pt > 1. && quality(\"highPurity\") && "
                        "hitPattern().numberOfLostHits(\"MISSING_OUTER_HITS\") < 5");
  desc.add<double>("delta_tk_ts_layer1", 0.02);
  desc.add<double>("delta_tk_ts_interface", 0.03);
  desc.add<double>("delta_ts_em_had", 0.03);
  desc.add<double>("delta_ts_had_had", 0.03);
  desc.add<double>("separationSmall", 2);  //cm
  desc.add<double>("separation", 6);       //cm
  desc.add<int>("maxDepth", 10);
  desc.add<double>("track_time_quality_threshold", 0.5);
  LinkingAlgoBase::fillPSetDescription(desc);
}
