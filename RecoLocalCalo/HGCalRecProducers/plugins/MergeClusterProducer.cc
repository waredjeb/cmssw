// Authors: Olivie Franklova - olivie.abigail.franklova@cern.ch
// Date: 03/2023
// @file merge layer clusters which were produce by HGCalLayerClusterProducer

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimCalorimetry/HGCalAssociatorProducers/interface/AssociatorTools.h"
#include "PreliminaryClusterFilterBase.h"
#include "PreliminaryClusterFilterFactory.h"
#include "PreliminaryClusterFilterByCP.h"

class MergeClusterProducer : public edm::stream::EDProducer<> {
public:
  /**
   * @brief Constructor with parameter settings - which can be changed in  ...todo.
   * Constructor will set all variables by input param ps. 
   * 
   * @param[in] ps parametr set to set variables
  */
  MergeClusterProducer(const edm::ParameterSet &);
  ~MergeClusterProducer() override {}
  /**
   * @brief Method fill description which will be used in pyhton file.
   * 
   * @param[out] description to be fill
  */
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  /**
   * @brief Method will merge the producers and put them back to event
   * 
   * @param[in, out] evt from get info and put result
   * @param[in] es to get event setup info
  */
  void produce(edm::Event &, const edm::EventSetup &) override;

private:
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> EEclusters_token_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> HSiclusters_token_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> HSciclusters_token_;
  edm::EDGetTokenT<std::vector<CaloParticle>> caloparticles_token_;

  std::string timeClname_;
  std::string clusterFilter_;
  const edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> clustersTimeEE_token_;
  const edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> clustersTimeHSi_token_;
  const edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> clustersTimeHSci_token_;
  std::unique_ptr<const ticl::PreliminaryClusterFilterBase> theFilter_;

  /**
   * @brief method merge three vectors of reco::CaloCluster to one
   * 
   * @param[out] merge the vector into which others vectors will be merge
   * @param[in] EE vector for Electromagnetic silicon
   * @param[in] HSi vector for Hardon silicon
   * @param[in] ESci vector for hadron scintillator
  */
  void mergeTogether(std::vector<reco::CaloCluster> &merge,
                     const std::vector<reco::CaloCluster> &EE,
                     const std::vector<reco::CaloCluster> &HSi,
                     const std::vector<reco::CaloCluster> &HSci);

  /**
   * @brief copy all values from vm to to
   * 
   * @param[in] vm Value map with values
   * @param[out] to vector to will be copy value map
  */
  void addTo(std::vector<std::pair<float, float>> &to,
             const edm::ValueMap<std::pair<float, float>> &vm,
             const std::vector<float> &layerClusterMask) {
    for (size_t i = 0; i < layerClusterMask.size(); ++i) {
      if (layerClusterMask[i] == 1.f) {
        to.push_back(vm.get(i));
      }
    }
  }
  /**
   * @brief Merge value map of time for all parts of detector together  to vector times
   * 
   * @param[in] evt Event to get time value maps
   * @param[in] size of all 3 value maps
   * @param[out] times vector of merged time vectors
  */
  void mergeTime(edm::Event &evt,
                 size_t size,
                 std::vector<std::pair<float, float>> &times,
                 const std::vector<float> &layerClusterMask) {
    edm::Handle<edm::ValueMap<std::pair<float, float>>> EE, HSi, HSci;
    // get values from all three part of detectors
    evt.getByToken(clustersTimeEE_token_, EE);
    evt.getByToken(clustersTimeHSi_token_, HSi);
    evt.getByToken(clustersTimeHSci_token_, HSci);
    addTo(times, *EE, layerClusterMask);
    addTo(times, *HSi, layerClusterMask);
    addTo(times, *HSci, layerClusterMask);
  }
  /**
   * @brief get info form event and then call merge
   * 
   * it is used for merge and clusters and time
   * 
   * @param[in] evt Event
   * @param[in] EE_token token for Electromagnetic silicon
   * @param[in] HSi_token token for Hardon silicon
   * @param[in] ESci_token token for hadron scintillator
   * @return merged result
  */
  template <typename T>
  void createMerge(edm::Event &evt,
                   const edm::EDGetTokenT<T> &EE_token,
                   const edm::EDGetTokenT<T> &HSi_token,
                   const edm::EDGetTokenT<T> &HSci_token,
                   T &merge) {
    edm::Handle<T> EE, HSi, HSci;
    // get values from all three part of detectors
    evt.getByToken(EE_token, EE);
    evt.getByToken(HSi_token, HSi);
    evt.getByToken(HSci_token, HSci);
    mergeTogether(merge, *EE, *HSi, *HSci);
  }
};

MergeClusterProducer::MergeClusterProducer(const edm::ParameterSet &ps)
    : timeClname_(ps.getParameter<std::string>("timeClname")),
      clustersTimeEE_token_(
          consumes<edm::ValueMap<std::pair<float, float>>>(ps.getParameter<edm::InputTag>("time_layerclustersEE"))),
      clustersTimeHSi_token_(
          consumes<edm::ValueMap<std::pair<float, float>>>(ps.getParameter<edm::InputTag>("time_layerclustersHSi"))),
      clustersTimeHSci_token_(
          consumes<edm::ValueMap<std::pair<float, float>>>(ps.getParameter<edm::InputTag>("time_layerclustersHSci"))) {
  EEclusters_token_ = consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layerClustersEE"));
  HSiclusters_token_ = consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layerClustersHSi"));
  HSciclusters_token_ = consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layerClustersHSci"));
  caloparticles_token_ = consumes<std::vector<CaloParticle>>(ps.getParameter<edm::InputTag>("caloparticles"));
  clusterFilter_ = ps.getParameter<std::string>("clusterFilter");
  theFilter_ = PreliminaryClusterFilterFactory::get()->create(clusterFilter_, ps);

  produces<std::vector<float>>("InitialLayerClustersMask");
  produces<std::vector<reco::BasicCluster>>();
  produces<std::vector<reco::BasicCluster>>("sharing");
  //time for layer clusters
  produces<edm::ValueMap<std::pair<float, float>>>(timeClname_);
}

void MergeClusterProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // hgcalMergeLayerClusters
  edm::ParameterSetDescription desc;
  //layer clusters
  desc.add<edm::InputTag>("layerClustersEE", edm::InputTag("hgcalLayerClustersEE"));
  desc.add<edm::InputTag>("layerClustersHSi", edm::InputTag("hgcalLayerClustersHSi"));
  desc.add<edm::InputTag>("layerClustersHSci", edm::InputTag("hgcalLayerClustersHSci"));
  desc.add<edm::InputTag>("caloparticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<std::string>("clusterFilter", "PreliminaryClusterFilterByCP");
  desc.add<double>("deltaEta", 0.3);
  desc.add<double>("deltaPhi", 1.57);
  //time
  desc.add<edm::InputTag>("time_layerclustersEE", edm::InputTag("hgcalLayerClustersEE", "timeLayerCluster"));
  desc.add<edm::InputTag>("time_layerclustersHSi", edm::InputTag("hgcalLayerClustersHSi", "timeLayerCluster"));
  desc.add<edm::InputTag>("time_layerclustersHSci", edm::InputTag("hgcalLayerClustersHSci", "timeLayerCluster"));

  desc.add<std::string>("timeClname", "timeLayerCluster");
  descriptions.add("hgcalMergeLayerClusters", desc);
}

void MergeClusterProducer::produce(edm::Event &evt, const edm::EventSetup &es) {
  //merge clusters
  std::unique_ptr<std::vector<reco::BasicCluster>> clusters(new std::vector<reco::BasicCluster>);
  auto clustersFiltered = std::make_unique<std::vector<reco::BasicCluster>>();
  auto const &caloParticles = evt.get(caloparticles_token_);
  std::vector<size_t> caloParticleIndices;
  removeCPFromPU(caloParticles, caloParticleIndices, true);

  createMerge(evt, EEclusters_token_, HSiclusters_token_, HSciclusters_token_, *clusters);

  std::unique_ptr<std::vector<float>> layerClustersMask(new std::vector<float>(clusters->size(), 0.f));
  theFilter_->filter(*clusters, *layerClustersMask, caloParticles, caloParticleIndices);
  std::unique_ptr<std::vector<float>> layerClustersMaskFiltered(new std::vector<float>);
  edm::Handle<edm::ValueMap<std::pair<float, float>>> EE, HSi, HSci;
  // get values from all three part of detectors
  evt.getByToken(clustersTimeEE_token_, EE);
  evt.getByToken(clustersTimeHSi_token_, HSi);
  evt.getByToken(clustersTimeHSci_token_, HSci);

  edm::Handle<std::vector<reco::BasicCluster>> EEC, HSiC, HSciC;
  // get values from all three part of detectors
  evt.getByToken(EEclusters_token_, EEC);
  evt.getByToken(HSiclusters_token_, HSiC);
  evt.getByToken(HSciclusters_token_, HSciC);

  std::vector<std::pair<float, float>> times;

  for (size_t i = 0; i < layerClustersMask->size(); i++) {
    if ((*layerClustersMask)[i] == 1.f) {
      clustersFiltered->push_back((*clusters)[i]);
      layerClustersMaskFiltered->push_back(1.f);
      if (std::find(EEC->begin(), EEC->end(), (*clusters)[i]) != EEC->end()) {
        times.push_back(EE->get(i));
      } else if (std::find(HSiC->begin(), HSiC->end(), (*clusters)[i]) != HSiC->end()) {
        times.push_back(HSi->get(i));
      } else if (std::find(HSciC->begin(), HSciC->end(), (*clusters)[i]) != HSciC->end()) {
        times.push_back(HSci->get(i));
      }
    }
  }
  //put new clusters to event
  auto clusterHandle = evt.put(std::move(clustersFiltered));

  //create layer cluster mask

  //put it into event
  evt.put(std::move(layerClustersMaskFiltered), "InitialLayerClustersMask");

  //time
  //  mergeTime(evt, clusterHandle->size(), times, *layerClustersMask);
  auto timeCl = std::make_unique<edm::ValueMap<std::pair<float, float>>>();
  edm::ValueMap<std::pair<float, float>>::Filler filler(*timeCl);
  filler.insert(clusterHandle, times.begin(), times.end());
  filler.fill();
  evt.put(std::move(timeCl), timeClname_);
}

void MergeClusterProducer::mergeTogether(std::vector<reco::CaloCluster> &merge,
                                         const std::vector<reco::CaloCluster> &EE,
                                         const std::vector<reco::CaloCluster> &HSi,
                                         const std::vector<reco::CaloCluster> &HSci) {
  auto clusterSize = EE.size() + HSi.size() + HSci.size();
  merge.reserve(clusterSize);

  merge.insert(merge.end(), EE.begin(), EE.end());
  merge.insert(merge.end(), HSi.begin(), HSi.end());
  merge.insert(merge.end(), HSci.begin(), HSci.end());
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MergeClusterProducer);
