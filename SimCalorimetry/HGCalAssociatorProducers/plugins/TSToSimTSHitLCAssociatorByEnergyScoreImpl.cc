#include <cfloat>

#include "TSToSimTSHitLCAssociatorByEnergyScoreImpl.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimCalorimetry/HGCalAssociatorProducers/interface/AssociatorTools.h"

TSToSimTSHitLCAssociatorByEnergyScoreImpl::TSToSimTSHitLCAssociatorByEnergyScoreImpl(
    edm::EDProductGetter const& productGetter,
    bool hardScatterOnly,
    std::shared_ptr<hgcal::RecHitTools> recHitTools,
    const std::unordered_map<DetId, const HGCRecHit*>* hitMap)
    : hardScatterOnly_(hardScatterOnly), recHitTools_(recHitTools), hitMap_(hitMap), productGetter_(&productGetter) {
  layers_ = recHitTools_->lastLayerBH();
}

hgcal::association TSToSimTSHitLCAssociatorByEnergyScoreImpl::makeConnections(
    const edm::Handle<ticl::TracksterCollection>& tCH,
    const edm::Handle<reco::CaloClusterCollection>& lCCH,
    const edm::Handle<SimClusterCollection>& sCCH,
    const edm::Handle<CaloParticleCollection>& cPCH,    
    const edm::Handle<std::map<uint, std::vector<uint>>>& simTrackstersMapH,
    const edm::Handle<ticl::TracksterCollection>& sTCH,
    const edm::Handle<ticl::TracksterCollection>& sTfromCPCH,
    const hgcal::validationType valType) const {
  // Get collections
  const auto& tracksters = *tCH.product();
  const auto& layerClusters = *lCCH.product();
  const auto& sC = *sCCH.product();                              
  const auto& cP = *cPCH.product();
  const auto cPHandle_id = cPCH.id();
  const auto& cpToSc_SimTrackstersMap = *simTrackstersMapH.product();
  const auto& simTSs = *sTCH.product();
  const auto& simTSs_fromCP = *sTfromCPCH.product();
  const auto nSC = sC.size();
  const auto nTracksters = tracksters.size();
  const auto nSimTracksters = simTSs.size();

  std::unordered_map<DetId, std::vector<hgcal::detIdInfoInCluster>> detIdSimTSId_Map;
  std::unordered_map<DetId, std::vector<hgcal::detIdInfoInTrackster>> detIdToTracksterId_Map;

  // This vector contains the ids of the SimTracksters contributing with at least one hit to the Trackster and the association score
  //stsInTrackster[trackster][STSids]
  //Connects a Trackster with all related SimTracksters.
  hgcal::tracksterToSimTrackster stsInTrackster;  // tsId->(stId,score)
  stsInTrackster.resize(nTracksters);

  hgcal::simTracksterToTrackster tssInSimTrackster;
  tssInSimTrackster.resize(nSimTracksters);

  // cPOnLayer[caloparticle]:
  //1. the sum of all rechits energy times fraction of the relevant simhit related to that calo particle.
  //2. the hits and fractions of that calo particle.
  //3. the layer clusters with matched rechit id.
  std::unordered_map<int, hgcal::caloParticleOnLayer> cPOnLayer;
  std::unordered_map<int, std::vector<hgcal::caloParticleOnLayer>> sCOnLayer;
  //Consider CaloParticles coming from the hard scatterer, excluding the PU contribution.
  std::vector<size_t> cPIndices;
  removeCPFromPU(cP, cPIndices, hardScatterOnly_);
  for (const auto cpIndex : cPIndices) {
    cPOnLayer[cpIndex].caloParticleId = cpIndex;
    cPOnLayer[cpIndex].energy = 0.f;
    cPOnLayer[cpIndex].hits_and_fractions.clear();
    sCOnLayer[cpIndex].resize(nSC);
    for (unsigned int iSC = 0; iSC < nSC; iSC++) {
      sCOnLayer[cpIndex][iSC].caloParticleId = cpIndex;
      sCOnLayer[cpIndex][iSC].energy = 0.f;
      sCOnLayer[cpIndex][iSC].hits_and_fractions.clear();
    }
  }

  auto getCPId = [](const ticl::Trackster& simTS,
                    const unsigned int iSTS,
                    const edm::ProductID& cPHandle_id,
                    const std::map<unsigned int, std::vector<unsigned int>>& cpToSc_SimTrackstersMap,
                    const ticl::TracksterCollection& simTSs_fromCP) {
    unsigned int cpId = -1;

    const auto productID = simTS.seedID();
    if (productID == cPHandle_id) {
      cpId = simTS.seedIndex();
    } else {  // SimTrackster from SimCluster
      const auto findSimTSFromCP = std::find_if(
          std::begin(cpToSc_SimTrackstersMap),
          std::end(cpToSc_SimTrackstersMap),
          [&](const std::pair<unsigned int, std::vector<unsigned int>>& cpToScs) {
            return std::find(std::begin(cpToScs.second), std::end(cpToScs.second), iSTS) != std::end(cpToScs.second);
          });
      if (findSimTSFromCP != std::end(cpToSc_SimTrackstersMap)) {
        cpId = simTSs_fromCP[findSimTSFromCP->first].seedIndex();
      }
    }

    return cpId;
  };

  auto getLCId = [](const std::vector<unsigned int>& tst_vertices,
                    const reco::CaloClusterCollection& layerClusters,
                    const DetId& hitid) {
    unsigned int lcId = -1;
    std::for_each(std::begin(tst_vertices), std::end(tst_vertices), [&](unsigned int idx) {
      const auto& lc_haf = layerClusters[idx].hitsAndFractions();
      const auto& hitFound = std::find_if(std::begin(lc_haf),
                                          std::end(lc_haf),
                                          [&hitid](const std::pair<DetId, float>& v) { return v.first == hitid; });
      if (hitFound != lc_haf.end())  // not all hits may be clusterized
        lcId = idx;
    });
    return lcId;
  };

  for (unsigned int iSTS = 0; iSTS < nSimTracksters; ++iSTS) {
    const auto cpId = getCPId(simTSs[iSTS], iSTS, cPHandle_id, cpToSc_SimTrackstersMap, simTSs_fromCP);
    if (std::find(cPIndices.begin(), cPIndices.end(), cpId) == cPIndices.end())
      continue;

    // Loop through SimClusters
    for (const auto& simCluster : cP[cpId].simClusters()) {
      auto iSim = simTSs[iSTS].seedIndex();
      if (simTSs[iSTS].seedID() != cPHandle_id) {  // SimTrackster from SimCluster
        if (iSim != (&(*simCluster) - &(sC[0])))
          continue;
      } else
        iSim = 0;

      for (const auto& it_haf : simCluster->hits_and_fractions()) {
        const auto hitid = (it_haf.first);
        const auto lcId = getLCId(simTSs[iSTS].vertices(), layerClusters, hitid);
        //V9:maps the layers in -z: 0->51 and in +z: 52->103
        //V10:maps the layers in -z: 0->49 and in +z: 50->99
        const auto itcheck = hitMap_->find(hitid);
        //Checks whether the current hit belonging to sim cluster has a reconstructed hit.
        if ((valType == 0 && itcheck != hitMap_->end()) || (valType > 0 && int(lcId) >= 0)) {
          const auto elemId = (valType == 0) ? hitid : lcId;
          const auto iLC = std::find(simTSs[iSTS].vertices().begin(), simTSs[iSTS].vertices().end(), lcId);
          const auto lcFraction =
              1.f / simTSs[iSTS].vertex_multiplicity(std::distance(std::begin(simTSs[iSTS].vertices()), iLC));
          const auto elemFr = (valType == 0) ? it_haf.second : lcFraction;
          //Since the current hit from sim cluster has a reconstructed hit with the same detid,
          //make a map that will connect a detid with:
          //1. the CaloParticles that have a SimCluster with sim hits in that cell via caloparticle id.
          //2. the sum of all SimHits fractions that contributes to that detid.
          //So, keep in mind that in case of multiple CaloParticles contributing in the same cell
          //the fraction is the sum over all calo particles. So, something like:
          //detid: (caloparticle 1, sum of hits fractions in that detid over all cp) , (caloparticle 2, sum of hits fractions in that detid over all cp), (caloparticle 3, sum of hits fractions in that detid over all cp) ...
          if (detIdSimTSId_Map.find(elemId) == detIdSimTSId_Map.end()) {
            detIdSimTSId_Map[elemId] = std::vector<hgcal::detIdInfoInCluster>();
            detIdSimTSId_Map[elemId].emplace_back(hgcal::detIdInfoInCluster{iSTS, elemFr});
          } else {
            auto findSTSIt =
                std::find(detIdSimTSId_Map[elemId].begin(),
                          detIdSimTSId_Map[elemId].end(),
                          hgcal::detIdInfoInCluster{
                              iSTS, 0});  // only the first element is used for the matching (overloaded operator==)
            if (findSTSIt != detIdSimTSId_Map[elemId].end()) {
              if (valType == 0)
                findSTSIt->fraction += elemFr;
            } else {
              detIdSimTSId_Map[elemId].emplace_back(hgcal::detIdInfoInCluster{iSTS, elemFr});
            }
          }
          const auto hitEn = itcheck->second->energy();
          //Since the current hit from sim cluster has a reconstructed hit with the same detid,
          //fill the cPOnLayer[caloparticle][layer] object with energy (sum of all rechits energy times fraction
          //of the relevant simhit) and keep the hit (detid and fraction) that contributed.
          cPOnLayer[cpId].energy += it_haf.second * hitEn;
          sCOnLayer[cpId][iSim].energy += lcFraction * hitEn;
          // Need to compress the hits and fractions in order to have a
          // reasonable score between CP and LC. Imagine, for example, that a
          // CP has detID X used by 2 SimClusters with different fractions. If
          // a single LC uses X with fraction 1 and is compared to the 2
          // contributions separately, it will be assigned a score != 0, which
          // is wrong.
          auto& haf = cPOnLayer[cpId].hits_and_fractions;
          auto found = std::find_if(
              std::begin(haf), std::end(haf), [&hitid](const std::pair<DetId, float>& v) { return v.first == hitid; });
          if (found != haf.end())
            found->second += it_haf.second;
          else
            haf.emplace_back(hitid, it_haf.second);
          // Same for sCOnLayer
          auto& haf_sc = sCOnLayer[cpId][iSim].hits_and_fractions;
          auto found_sc = std::find_if(std::begin(haf_sc),
                                       std::end(haf_sc),
                                       [&hitid](const std::pair<DetId, float>& v) { return v.first == hitid; });
          if (found_sc != haf_sc.end())
            found_sc->second += it_haf.second;
          else
            haf_sc.emplace_back(hitid, it_haf.second);
        }
      }  // end of loop through SimHits
    }    // end of loop through SimClusters
  }      // end of loop through SimTracksters

  auto apply_LCMultiplicity = [](const ticl::Trackster& trackster, const reco::CaloClusterCollection& layerClusters) {
    std::vector<std::pair<DetId, float>> hits_and_fractions_norm;
    int lcInTst = 0;
    std::for_each(std::begin(trackster.vertices()), std::end(trackster.vertices()), [&](unsigned int idx) {
      const auto fraction = 1.f / trackster.vertex_multiplicity(lcInTst++);
      for (const auto& cell : layerClusters[idx].hitsAndFractions()) {
        hits_and_fractions_norm.emplace_back(
            cell.first, cell.second * fraction);  // cell.second is the hit fraction in the layerCluster
      }
    });
    return hits_and_fractions_norm;
  };

  // Loop through Tracksters
  std::unordered_map<int, std::pair<float, float>> tracksterIdToEnergyAndScoreMap;
  for (unsigned int tstId = 0; tstId < nTracksters; ++tstId) {
    const auto& tst = tracksters[tstId];
    if (tst.vertices().empty())
      continue;

    const auto tst_hitsAndFractions = apply_LCMultiplicity(tst, layerClusters);
    const auto numberOfHitsInTS = tst_hitsAndFractions.size();
    //Loop through the hits of the trackster under study
    for (unsigned int iHit = 0; iHit < numberOfHitsInTS; iHit++) {
      const auto rh_detid = tst_hitsAndFractions[iHit].first;
      const auto rhFraction = tst_hitsAndFractions[iHit].second;

      const auto lcId_r = getLCId(tst.vertices(), layerClusters, rh_detid);
      const auto iLC_r = std::find(tst.vertices().begin(), tst.vertices().end(), lcId_r);
      const auto lcFraction_r = 1.f / tst.vertex_multiplicity(std::distance(std::begin(tst.vertices()), iLC_r));

      //Make a map that will connect a detid (that belongs to a rechit of the layer cluster under study,
      //no need to save others) with:
      //1. the layer clusters that have rechits in that detid
      //2. the fraction of the rechit of each layer cluster that contributes to that detid.
      //So, something like:
      //detid: (layer cluster 1, hit fraction) , (layer cluster 2, hit fraction), (layer cluster 3, hit fraction) ...
      //here comparing with the calo particle map above the
      if (detIdToTracksterId_Map.find(rh_detid) == detIdToTracksterId_Map.end()) {
        detIdToTracksterId_Map[rh_detid] = std::vector<hgcal::detIdInfoInTrackster>();
        detIdToTracksterId_Map[rh_detid].emplace_back(
            hgcal::detIdInfoInTrackster{tstId, lcId_r, rhFraction});
      } else {
        auto findTSIt =
            std::find(detIdToTracksterId_Map[rh_detid].begin(),
                      detIdToTracksterId_Map[rh_detid].end(),
                      hgcal::detIdInfoInTrackster{
                          tstId, 0, 0});  // only the first element is used for the matching (overloaded operator==)
        if (findTSIt != detIdToTracksterId_Map[rh_detid].end()) {
          if (valType == 0)
            findTSIt->fraction += rhFraction;
        } else {
          detIdToTracksterId_Map[rh_detid].emplace_back(
              hgcal::detIdInfoInTrackster{tstId, lcId_r, rhFraction});
        }
      }

      // Check whether the RecHit of the trackster under study has a SimHit in the same cell
      const auto elemId = (valType == 0) ? rh_detid.rawId() : lcId_r;
      const auto recoFr = (valType == 0) ? rhFraction : lcFraction_r;
      const auto& hit_find_in_STS = detIdSimTSId_Map.find(elemId);
      if (hit_find_in_STS != detIdSimTSId_Map.end()) {
        // Since the hit is belonging to the layer cluster, it must be also in the rechits map
        const auto hitEn = hitMap_->find(rh_detid)->second->energy();

        for (const auto& h : hit_find_in_STS->second) {
          const auto shared_fraction = std::min(recoFr, h.fraction);
          const auto iSTS = h.clusterId;
          const auto& simTS = simTSs[iSTS];
          auto iSim = simTS.seedIndex();
          if (simTSs[iSTS].seedID() == cPHandle_id)  // SimTrackster from CaloParticle
            iSim = 0;

          // SimTrackster with simHits connected via detid with the rechit under study
          //So, from all layers clusters, find the rechits that are connected with a calo particle and save/calculate the
          //energy of that calo particle as the sum over all rechits of the rechits energy weighted
          //by the caloparticle's fraction related to that rechit.
          const auto cpId = getCPId(simTS, iSTS, cPHandle_id, cpToSc_SimTrackstersMap, simTSs_fromCP);
          if (std::find(cPIndices.begin(), cPIndices.end(), cpId) == cPIndices.end())
            continue;

          tssInSimTrackster[iSTS][tstId].first += shared_fraction * hitEn;
          //Here cPOnLayer[caloparticle][layer] describe above is set.
          //Here for Tracksters with matched rechit the CP fraction times hit energy is added and saved
          cPOnLayer[cpId].layerClusterIdToEnergyAndScore[tstId].first += shared_fraction * hitEn;
          sCOnLayer[cpId][iSim].layerClusterIdToEnergyAndScore[tstId].first += shared_fraction * hitEn;
          cPOnLayer[cpId].layerClusterIdToEnergyAndScore[tstId].second = FLT_MAX;
          sCOnLayer[cpId][iSim].layerClusterIdToEnergyAndScore[tstId].second = FLT_MAX;
          stsInTrackster[tstId].emplace_back(iSTS, FLT_MAX);
        }
        
      }

    }  //end of loop through rechits of the layer cluster.

  }  // end of loop through Tracksters

  // Loop through Tracksters
  for (unsigned int tstId = 0; tstId < nTracksters; ++tstId) {
    const auto& tst = tracksters[tstId];
    if (tst.vertices().empty())
      continue;

    // find the unique SimTrackster ids contributing to the Trackster
    std::sort(stsInTrackster[tstId].begin(), stsInTrackster[tstId].end());
    const auto last = std::unique(stsInTrackster[tstId].begin(), stsInTrackster[tstId].end());
    stsInTrackster[tstId].erase(last, stsInTrackster[tstId].end());

    if (tst.raw_energy() == 0. && !stsInTrackster[tstId].empty()) {
      continue;
    }

    const auto tst_hitsAndFractions = apply_LCMultiplicity(tst, layerClusters);

    // Compute the correct normalization
    float tracksterEnergy = 0.f, invTracksterEnergyWeight = 0.f;
    for (const auto& haf : tst_hitsAndFractions) {
      float hitFr = 0.f;
      if (valType == 0) {
        hitFr = haf.second;
      } else {
        const auto lcId = getLCId(tst.vertices(), layerClusters, haf.first);
        const auto iLC = std::find(tst.vertices().begin(), tst.vertices().end(), lcId);
        hitFr = 1.f / tst.vertex_multiplicity(std::distance(std::begin(tst.vertices()), iLC));
      }
      tracksterEnergy += hitFr * hitMap_->at(haf.first)->energy();
      invTracksterEnergyWeight += pow(hitFr * hitMap_->at(haf.first)->energy(), 2);
    }
    if (invTracksterEnergyWeight)
      invTracksterEnergyWeight = 1.f / invTracksterEnergyWeight;

    for (const auto& haf : tst_hitsAndFractions) {
      const auto rh_detid = haf.first;
      unsigned int elemId = 0;
      float rhFraction = 0.f;
      if (valType == 0) {
        elemId = rh_detid.rawId();
        rhFraction = haf.second;
      } else {
        const auto lcId = getLCId(tst.vertices(), layerClusters, rh_detid);
        elemId = lcId;
        const auto iLC = std::find(tst.vertices().begin(), tst.vertices().end(), lcId);
        rhFraction = 1.f / tst.vertex_multiplicity(std::distance(std::begin(tst.vertices()), iLC));
      }

      bool hitWithNoSTS = false;
      if (detIdSimTSId_Map.find(elemId) == detIdSimTSId_Map.end())
        hitWithNoSTS = true;
      const HGCRecHit* hit = hitMap_->find(rh_detid)->second;
      const auto hitEnergyWeight = pow(hit->energy(), 2);

      for (auto& stsPair : stsInTrackster[tstId]) {
        float cpFraction = 0.f;
        if (!hitWithNoSTS) {
          const auto& findSTSIt = std::find(
              detIdSimTSId_Map[elemId].begin(),
              detIdSimTSId_Map[elemId].end(),
              hgcal::detIdInfoInCluster{
                  stsPair.first, 0.f});  // only the first element is used for the matching (overloaded operator==)
          if (findSTSIt != detIdSimTSId_Map[elemId].end())
            cpFraction = findSTSIt->fraction;
        }
        if (stsPair.second == FLT_MAX) {
          stsPair.second = 0.f;
        }
        stsPair.second +=
            std::min(pow(rhFraction - cpFraction, 2), pow(rhFraction, 2)) * hitEnergyWeight * invTracksterEnergyWeight;
      }
    }  // end of loop through trackster rechits
  }  // end of loop through Tracksters

  std::unordered_map<unsigned int, std::vector<float>> score3d;
  std::unordered_map<unsigned int, std::vector<float>> tstSharedEnergy;

  for (unsigned int iSTS = 0; iSTS < nSimTracksters; ++iSTS) {
    score3d[iSTS].resize(nTracksters);
    tstSharedEnergy[iSTS].resize(nTracksters);
    for (unsigned int j = 0; j < nTracksters; ++j) {
      score3d[iSTS][j] = FLT_MAX;
      tstSharedEnergy[iSTS][j] = 0.f;
    }
  }

  for (unsigned int iSTS = 0; iSTS < nSimTracksters; ++iSTS) {
    const auto& sts = simTSs[iSTS];
    const auto& cpId = getCPId(sts, iSTS, cPHandle_id, cpToSc_SimTrackstersMap, simTSs_fromCP);
    if (valType == 0 && std::find(cPIndices.begin(), cPIndices.end(), cpId) == cPIndices.end())
      continue;

    const auto& hafLC = apply_LCMultiplicity(sts, layerClusters);
    float SimEnergy_LC = 0.f;
    for (const auto& haf : hafLC) {
      const auto lcId = getLCId(sts.vertices(), layerClusters, haf.first);
      const auto iLC = std::find(sts.vertices().begin(), sts.vertices().end(), lcId);
      SimEnergy_LC +=
          hitMap_->at(haf.first)->energy() / sts.vertex_multiplicity(std::distance(std::begin(sts.vertices()), iLC));
    }

    auto iSim = sts.seedIndex();
    if (sts.seedID() == cPHandle_id)  // SimTrackster from CaloParticle
      iSim = 0;
    auto& simOnLayer = (valType == 0) ? cPOnLayer[cpId] : sCOnLayer[cpId][iSim];

    // Keep the Trackster ids that are related to
    // SimTrackster under study for the final filling of the score
    std::set<unsigned int> stsId_tstId_related;
    auto& score3d_iSTS = score3d[iSTS];

    float SimEnergy = 0.f;
    float SimEnergyWeight = 0.f, hitsEnergyWeight = 0.f;
    //for (unsigned int layerId = 0; layerId < 1/*layers * 2*/; ++layerId) {
    const auto SimNumberOfHits = simOnLayer.hits_and_fractions.size();
    if (SimNumberOfHits == 0)
      continue;
    SimEnergy += simOnLayer.energy;
 
    for (const auto& haf : ((valType == 0) ? simOnLayer.hits_and_fractions : hafLC)) {
      const auto& hitDetId = haf.first;
      // Compute the correct normalization
      // Need to loop on the simOnLayer data structure since this is the
      // only one that has the compressed information for multiple usage
      // of the same DetId by different SimClusters by a single CaloParticle.
      SimEnergyWeight += pow(haf.second * hitMap_->at(hitDetId)->energy(), 2);

      const auto lcId = getLCId(sts.vertices(), layerClusters, hitDetId);
      float cpFraction = 0.f;
      if (valType == 0) {
        cpFraction = haf.second;
      } else {
        const auto iLC = std::find(sts.vertices().begin(), sts.vertices().end(), lcId);
        cpFraction = 1.f / sts.vertex_multiplicity(std::distance(std::begin(sts.vertices()), iLC));
      }
      if (cpFraction == 0.f)
        continue;  // hopefully this should never happen

      bool hitWithNoTS = false;
      if (detIdToTracksterId_Map.find(hitDetId) == detIdToTracksterId_Map.end())
        hitWithNoTS = true;
      const HGCRecHit* hit = hitMap_->find(hitDetId)->second;
      const auto hitEnergyWeight = pow(hit->energy(), 2);
      hitsEnergyWeight += pow(cpFraction, 2) * hitEnergyWeight;

      for (auto& tsPair : simOnLayer.layerClusterIdToEnergyAndScore) {
        const auto tstId = tsPair.first;
        stsId_tstId_related.insert(tstId);

        float tstFraction = 0.f;
        if (!hitWithNoTS) {
          const auto findTSIt =
              std::find(detIdToTracksterId_Map[hitDetId].begin(),
                        detIdToTracksterId_Map[hitDetId].end(),
                        hgcal::detIdInfoInTrackster{
                            tstId, 0, 0.f});  // only the first element is used for the matching (overloaded operator==)
          if (findTSIt != detIdToTracksterId_Map[hitDetId].end()) {
            if (valType == 0) {
              tstFraction = findTSIt->fraction;
            } else {
              const auto iLC = std::find(
                  tracksters[tstId].vertices().begin(), tracksters[tstId].vertices().end(), findTSIt->clusterId);
              if (iLC != tracksters[tstId].vertices().end()) {
                tstFraction = 1.f / tracksters[tstId].vertex_multiplicity(
                                        std::distance(std::begin(tracksters[tstId].vertices()), iLC));
              }
            }
          }
        }
        // Here do not divide as before by the trackster energy weight. Should sum first
        // over all layers and divide with the total CP energy over all layers.
        if (tsPair.second.second == FLT_MAX) {
          tsPair.second.second = 0.f;
        }
        tsPair.second.second += std::min(pow(tstFraction - cpFraction, 2), pow(cpFraction, 2)) * hitEnergyWeight;

        LogDebug("HGCalValidator") << "\nTracksterId:\t" << tstId << "\tSimTracksterId:\t" << iSTS << "\tcpId:\t"
                                   << cpId << "\ttstfraction, cpfraction:\t" << tstFraction << ", " << cpFraction
                                   << "\thitEnergyWeight:\t" << hitEnergyWeight << "\tadded delta:\t"
                                   << pow((tstFraction - cpFraction), 2) * hitEnergyWeight
                                   << "\tcurrent Sim-score numerator:\t" << tsPair.second.second
                                   << "\tshared Sim energy:\t" << tsPair.second.first << '\n';
      }
    }  // end of loop through SimCluster SimHits on current layer

    if (simOnLayer.layerClusterIdToEnergyAndScore.empty())
      LogDebug("HGCalValidator") << "CP Id:\t" << cpId << "\tTS id:\t-1"
                                 << " Sub score in \t -1\n";

    for (const auto& tsPair : simOnLayer.layerClusterIdToEnergyAndScore) {
      const auto tstId = tsPair.first;
      // 3D score here without the denominator at this point
      if (score3d_iSTS[tstId] == FLT_MAX) {
        score3d_iSTS[tstId] = 0.f;
      }
      score3d_iSTS[tstId] += tsPair.second.second;
      tstSharedEnergy[iSTS][tstId] += tsPair.second.first;
    }

    const auto scoreDenom = (valType == 0) ? SimEnergyWeight : hitsEnergyWeight;
    const auto energyDenom = (valType == 0) ? SimEnergy : SimEnergy_LC;

    for (const auto tstId : stsId_tstId_related) {
      tssInSimTrackster[iSTS][tstId].first = tstSharedEnergy[iSTS][tstId] / energyDenom;
      tssInSimTrackster[iSTS][tstId].second = score3d_iSTS[tstId] /= scoreDenom;
    }
  }  // end of loop through SimTracksters


  return {stsInTrackster, tssInSimTrackster};
}

hgcal::RecoToSimCollectionSimTracksters TSToSimTSHitLCAssociatorByEnergyScoreImpl::associateRecoToSim(
    const edm::Handle<ticl::TracksterCollection>& tCH,
    const edm::Handle<reco::CaloClusterCollection>& lCCH,
    const edm::Handle<SimClusterCollection>& sCCH,
    const edm::Handle<CaloParticleCollection>& cPCH,
    const edm::Handle<std::map<uint, std::vector<uint>>>& simTrackstersMapH,
    const edm::Handle<ticl::TracksterCollection>& sTCH,
    const edm::Handle<ticl::TracksterCollection>& sTfromCPCH,
    const hgcal::validationType valType) const {
  hgcal::RecoToSimCollectionSimTracksters returnValue(productGetter_);
  const auto& links = makeConnections(tCH, lCCH, sCCH, cPCH, simTrackstersMapH, sTCH, sTfromCPCH, valType);

  const auto& stsInTrackster = std::get<0>(links);
  for (size_t tsId = 0; tsId < stsInTrackster.size(); ++tsId) {
    for (auto& stPair : stsInTrackster[tsId]) {
      LogDebug("TSToSimTSHitLCAssociatorByEnergyScoreImpl") << "Trackster Id:\t" << tsId << "\tSimTrackster id:\t"
                                                       << stPair.first << "\tscore:\t" << stPair.second << "\n";
      // Fill AssociationMap
      returnValue.insert(edm::Ref<ticl::TracksterCollection>(tCH, tsId),  // Ref to TS
                         std::make_pair(edm::Ref<ticl::TracksterCollection>(sTCH, stPair.first),
                                        stPair.second)  // Pair <Ref to ST, score>
      );
    }
  }
  return returnValue;
}

hgcal::SimToRecoCollectionSimTracksters TSToSimTSHitLCAssociatorByEnergyScoreImpl::associateSimToReco(
    const edm::Handle<ticl::TracksterCollection>& tCH,
    const edm::Handle<reco::CaloClusterCollection>& lCCH,
    const edm::Handle<SimClusterCollection>& sCCH,
    const edm::Handle<CaloParticleCollection>& cPCH,
    const edm::Handle<std::map<uint, std::vector<uint>>>& simTrackstersMapH,
    const edm::Handle<ticl::TracksterCollection>& sTCH,
    const edm::Handle<ticl::TracksterCollection>& sTfromCPCH,
    const hgcal::validationType valType) const {
  hgcal::SimToRecoCollectionSimTracksters returnValue(productGetter_);
  const auto& links = makeConnections(tCH, lCCH, sCCH, cPCH, simTrackstersMapH, sTCH, sTfromCPCH, valType);


  const auto& tssInSimTrackster = std::get<1>(links);
  for (size_t stId = 0; stId < tssInSimTrackster.size(); ++stId) {
    for (auto& tsPair : tssInSimTrackster[stId]) {

    LogDebug("TSToSimTSHitLCAssociatorByEnergyScoreImpl") << "SimTrackster Id:\t" << stId << "\tTrackster id:\t"
                                                       << tsPair.first << "\tscore:\t" << tsPair.second.second << "\n";
      returnValue.insert(
          edm::Ref<ticl::TracksterCollection>(sTCH, stId),                           // Ref to ST
          std::make_pair(edm::Ref<ticl::TracksterCollection>(tCH, tsPair.first),     // Pair <Ref to TS,
                         std::make_pair(tsPair.second.first, tsPair.second.second))  // pair <energy, score> >
      );
    }
  }
  return returnValue;
}
