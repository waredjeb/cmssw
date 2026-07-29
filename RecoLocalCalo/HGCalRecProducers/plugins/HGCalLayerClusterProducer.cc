// Authors: Olivie Franklova - olivie.abigail.franklova@cern.ch
// Date: 03/2023
// @file create layer clusters

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"

#include "RecoParticleFlow/PFClusterProducer/interface/RecHitTopologicalCleanerBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/SeedFinderBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/InitialClusteringStepBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterBuilderBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFCPositionCalculatorBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterEnergyCorrectorBase.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/ComputeClusterTime.h"

#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalLayerClusterAlgoFactory.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterHostCollection.h"
#include "DataFormats/TICL/interface/HitsAndFractionsHost.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"

#include <span>

class HGCalLayerClusterProducer : public edm::stream::EDProducer<> {
public:
  /**
   * @brief Constructor with parameter settings - which can be changed in hgcalLayerCluster_cff.py.
   * Constructor will set all variables by input param ps.
   * algoID variables will be set accordingly to the detector type.
   *
   * @param[in] ps parametr set to set variables
  */
  HGCalLayerClusterProducer(const edm::ParameterSet&);
  ~HGCalLayerClusterProducer() override {}
  /**
   * @brief Method fill description which will be used in pyhton file.
   *
   * @param[out] description to be fill
  */
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  /**
   * @brief Method run the algoritm to get clusters.
   *
   * @param[in, out] evt from get info and put result
   * @param[in] es to get event setup info
  */
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  edm::EDGetTokenT<HGCRecHitCollection> hits_token_;

  reco::CaloCluster::AlgoId algoId_;

  std::unique_ptr<HGCalClusteringAlgoBase> algo_;
  std::string detector_;

  unsigned int hitsTime_;

  // for calculate position
  std::vector<double> thresholdW0_;
  double positionDeltaRho2_;
  hgcal::RecHitTools rhtools_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  const bool calculatePositionInAlgo_;

  /**
   * @brief Sets algoId accordingly to the detector type
  */
  void setAlgoId();

  /**
   * @brief Counts position for all points in the cluster
   *
   * @param[in] hitmap hitmap to find correct RecHit
   * @param[in] hitsAndFractions hits of this cluster, from the association map
   * @return counted position
  */
  math::XYZPoint calculatePosition(std::unordered_map<uint32_t, const HGCRecHit*>& hitmap,
                                   std::span<const ticl::HitAndFraction> hitsAndFractions);

  /**
   * @brief Counts time for all points in the cluster
   *
   * @param[in] hitmap hitmap to find correct RecHit only for silicon (not for BH-HSci)
   * @param[in] hitsAndFractions hits of this cluster, from the association map
   * @return counted time
  */
  std::pair<float, float> calculateTime(std::unordered_map<uint32_t, const HGCRecHit*>& hitmap,
                                        std::span<const ticl::HitAndFraction> hitsAndFractions,
                                        size_t sizeCluster);
};

HGCalLayerClusterProducer::HGCalLayerClusterProducer(const edm::ParameterSet& ps)
    : algoId_(reco::CaloCluster::undefined),
      detector_(ps.getParameter<std::string>("detector")),  // one of EE, FH, BH, HFNose
      hitsTime_(ps.getParameter<unsigned int>("nHitsTime")),
      caloGeomToken_(consumesCollector().esConsumes<CaloGeometry, CaloGeometryRecord>()),
      calculatePositionInAlgo_(ps.getParameter<bool>("calculatePositionInAlgo")) {
  setAlgoId();  //sets algo id according to detector type
  hits_token_ = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("recHits"));

  auto pluginPSet = ps.getParameter<edm::ParameterSet>("plugin");
  if (detector_ == "HFNose") {
    algo_ = HGCalLayerClusterAlgoFactory::get()->create("HFNoseCLUE", pluginPSet);
    algo_->setAlgoId(algoId_, true);
  } else {
    algo_ = HGCalLayerClusterAlgoFactory::get()->create(pluginPSet.getParameter<std::string>("type"), pluginPSet);
    algo_->setAlgoId(algoId_);
  }
  thresholdW0_ = pluginPSet.getParameter<std::vector<double>>("thresholdW0");
  positionDeltaRho2_ = pluginPSet.getParameter<double>("positionDeltaRho2");

  produces<std::vector<float>>("InitialLayerClustersMask");

  produces<reco::CaloClusterHostCollection>();
  produces<ticl::HitsAndFractionsHost>();
}

void HGCalLayerClusterProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // hgcalLayerClusters
  edm::ParameterSetDescription desc;
  edm::ParameterSetDescription pluginDesc;
  pluginDesc.addNode(edm::PluginDescription<HGCalLayerClusterAlgoFactory>("type", "SiCLUE", true));

  desc.add<edm::ParameterSetDescription>("plugin", pluginDesc);
  desc.add<std::string>("detector", "EE")->setComment("options EE, FH, BH,  HFNose; other value defaults to EE");
  desc.add<edm::InputTag>("recHits", edm::InputTag("HGCalRecHit", "HGCEERecHits"));
  desc.add<std::string>("timeClname", "timeLayerCluster");
  desc.add<unsigned int>("nHitsTime", 3);
  desc.add<bool>("calculatePositionInAlgo", true);
  descriptions.add("hgcalLayerClusters", desc);
}

math::XYZPoint HGCalLayerClusterProducer::calculatePosition(
    std::unordered_map<uint32_t, const HGCRecHit*>& hitmap,
    std::span<const ticl::HitAndFraction> hitsAndFractions) {
  float total_weight = 0.f;
  float maxEnergyValue = 0.f;
  DetId maxEnergyIndex;
  float x = 0.f;
  float y = 0.f;

  for (auto const& hit : hitsAndFractions) {
    //time is computed wrt  0-25ns + offset and set to -1 if no time
    const HGCRecHit* rechit = hitmap[hit.hit];
    total_weight += rechit->energy();
    if (rechit->energy() > maxEnergyValue) {
      maxEnergyValue = rechit->energy();
      maxEnergyIndex = rechit->detid();
    }
  }
  float total_weight_log = 0.f;
  auto thick = rhtools_.getSiThickIndex(maxEnergyIndex);
  const GlobalPoint positionMaxEnergy(rhtools_.getPosition(maxEnergyIndex));
  for (auto const& hit : hitsAndFractions) {
    //time is computed wrt  0-25ns + offset and set to -1 if no time
    const HGCRecHit* rechit = hitmap[hit.hit];

    const GlobalPoint position(rhtools_.getPosition(rechit->detid()));

    if (thick != -1) {  //silicon
      //for silicon only just use 1+6 cells = 1.3cm for all thicknesses
      const float d1 = position.x() - positionMaxEnergy.x();
      const float d2 = position.y() - positionMaxEnergy.y();
      if ((d1 * d1 + d2 * d2) > positionDeltaRho2_)
        continue;

      float Wi = std::max(thresholdW0_[thick] + std::log(rechit->energy() / total_weight), 0.);
      x += position.x() * Wi;
      y += position.y() * Wi;
      total_weight_log += Wi;
    } else {  //scintillator
      x += position.x() * rechit->energy();
      y += position.y() * rechit->energy();
    }
  }
  if (thick != -1) {
    total_weight = total_weight_log;
  }
  if (total_weight != 0.) {
    float inv_tot_weight = 1.f / total_weight;
    return math::XYZPoint(x * inv_tot_weight, y * inv_tot_weight, positionMaxEnergy.z());
  } else {
    return {positionMaxEnergy.x(), positionMaxEnergy.y(), positionMaxEnergy.z()};
  }
}

std::pair<float, float> HGCalLayerClusterProducer::calculateTime(
    std::unordered_map<uint32_t, const HGCRecHit*>& hitmap,
    std::span<const ticl::HitAndFraction> hitsAndFractions,
    size_t sizeCluster) {
  std::pair<float, float> timeCl(-99., -1.);

  if (sizeCluster >= hitsTime_) {
    std::vector<float> timeClhits;
    std::vector<float> timeErrorClhits;

    for (auto const& hit : hitsAndFractions) {
      //time is computed wrt  0-25ns + offset and set to -1 if no time
      const HGCRecHit* rechit = hitmap[hit.hit];

      float rhTimeE = rechit->timeError();
      //check on timeError to exclude scintillator
      if (rhTimeE < 0.f)
        continue;
      timeClhits.push_back(rechit->time());
      timeErrorClhits.push_back(1.f / (rhTimeE * rhTimeE));
    }
    hgcalsimclustertime::ComputeClusterTime timeEstimator;
    timeCl = timeEstimator.fixSizeHighestDensity(timeClhits, timeErrorClhits, hitsTime_);
  }
  return timeCl;
}
void HGCalLayerClusterProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
  edm::Handle<HGCRecHitCollection> hits;

  edm::ESHandle<CaloGeometry> geom = es.getHandle(caloGeomToken_);
  rhtools_.setGeometry(*geom);
  algo_->getEventSetup(es, rhtools_);

  //make a map detid-rechit
  // NB for the moment just host EE and FH hits
  // timing in digi for BH not implemented for now
  std::unordered_map<uint32_t, const HGCRecHit*> hitmap;

  evt.getByToken(hits_token_, hits);
  algo_->populate(*hits);
  for (auto const& it : *hits) {
    hitmap[it.detid().rawId()] = &(it);
  }

  algo_->makeClusters();

  auto clustersAndAssociations = algo_->getClusters(false);
  auto clusters = std::move(clustersAndAssociations.layer_clusters);
  auto hitsAndFractions = std::move(clustersAndAssociations.hits_and_fractions);

  auto clusters_v = clusters->view();
  auto const hits_v = hitsAndFractions->const_view();
  const int numberOfClusters = clusters_v.position().metadata().size();

  for (int i = 0; i < numberOfClusters; ++i) {
    auto const hitsOfCluster = hits_v[i];
    if (!calculatePositionInAlgo_) {
      const auto position = calculatePosition(hitmap, hitsOfCluster);
      clusters_v.position()[i].x() = position.x();
      clusters_v.position()[i].y() = position.y();
      clusters_v.position()[i].z() = position.z();
    }
    // BH has no timing
    const auto timeCl = (detector_ != "BH") ? calculateTime(hitmap, hitsOfCluster, clusters_v.position()[i].cells())
                                            : std::pair<float, float>(-99.f, -1.f);
    clusters_v.timing()[i].time() = timeCl.first;
    clusters_v.timing()[i].timeError() = timeCl.second;
  }

  if (detector_ == "HFNose") {
    std::unique_ptr<std::vector<float>> layerClustersMask(new std::vector<float>);
    layerClustersMask->resize(numberOfClusters, 1.0);
    evt.put(std::move(layerClustersMask), "InitialLayerClustersMask");
  }

  evt.put(std::move(clusters));
  evt.put(std::move(hitsAndFractions));

  algo_->reset();
}

void HGCalLayerClusterProducer::setAlgoId() {
  if (detector_ == "EE") {
    algoId_ = reco::CaloCluster::hgcal_em;
  } else if (detector_ == "FH") {
    algoId_ = reco::CaloCluster::hgcal_had;
  } else if (detector_ == "BH") {
    algoId_ = reco::CaloCluster::hgcal_scintillator;
  } else if (detector_ == "HFNose") {
    algoId_ = reco::CaloCluster::hfnose;
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HGCalLayerClusterProducer);
