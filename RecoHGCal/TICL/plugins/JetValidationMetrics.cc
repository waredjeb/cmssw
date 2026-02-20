// JetValidationMetrics.cc

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"

#include <cmath>
#include <string>

class JetValidationMetrics : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit JetValidationMetrics(const edm::ParameterSet& iConfig);
  ~JetValidationMetrics() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  void endJob() override;

private:
  // inputs
  const edm::InputTag jetsTag_;
  const edm::InputTag genJetsTag_;

  const double recoJetPtThreshold_;
  const double matchGenPtThreshold_;
  const double rThreshold_;
  const double absEtaMax_;

  edm::EDGetTokenT<edm::View<reco::Jet>> jetsToken_;
  edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;

  // output
  TTree* outTree_ = nullptr;

  // accumulators 
  uint64_t nEvents_ = 0;
  uint64_t nGen_ = 0;
  uint64_t nMatched_ = 0;

  double sumResp_ = 0.0;
  double sumResp2_ = 0.0;

  double sumDR_ = 0.0;
  double sumDR2_ = 0.0;

   // low/high accumulators
   uint64_t nMatchedLowPt_ = 0;
   uint64_t nMatchedHighPt_ = 0;
 
   double sumRespLowPt_ = 0.0;
   double sumResp2LowPt_ = 0.0;
 
   double sumRespHighPt_ = 0.0;
   double sumResp2HighPt_ = 0.0;

  // tree branches 
  uint64_t b_nEvents_ = 0;
  uint64_t b_nGen_ = 0;
  uint64_t b_nMatched_ = 0;

  double b_meanRespLowPt_ = 0.0;
  double b_rmsRespLowPt_ = 0.0;
  double b_meanRespHighPt_ = 0.0;
  double b_rmsRespHighPt_ = 0.0;

};

JetValidationMetrics::JetValidationMetrics(const edm::ParameterSet& iConfig)
    : jetsTag_(iConfig.getParameter<edm::InputTag>("jets")),
      genJetsTag_(iConfig.getParameter<edm::InputTag>("genjets")),
      recoJetPtThreshold_(iConfig.getParameter<double>("recoJetPtThreshold")),
      matchGenPtThreshold_(iConfig.getParameter<double>("matchGenPtThreshold")),
      rThreshold_(iConfig.getParameter<double>("RThreshold")),
      absEtaMax_(iConfig.getParameter<double>("absEtaMax")) {
  jetsToken_ = consumes<edm::View<reco::Jet>>(jetsTag_);
  genJetsToken_ = consumes<reco::GenJetCollection>(genJetsTag_);


  edm::Service<TFileService> fs;
  const std::string label = moduleDescription().moduleLabel();
  TFileDirectory dir = fs->mkdir(label);

  outTree_ = dir.make<TTree>("output", "output");

  outTree_->Branch("nEvents", &b_nEvents_, "nEvents/l");
  outTree_->Branch("nGen", &b_nGen_, "nGen/l");
  outTree_->Branch("nMatched", &b_nMatched_, "nMatched/l");

  outTree_->Branch("meanRespLowPt", &b_meanRespLowPt_, "meanRespLowPt/D");
  outTree_->Branch("rmsRespLowPt", &b_rmsRespLowPt_, "rmsRespLowPt/D");
  outTree_->Branch("meanRespHighPt", &b_meanRespHighPt_, "meanRespHighPt/D");
  outTree_->Branch("rmsRespHighPt", &b_rmsRespHighPt_, "rmsRespHighPt/D");

}

void JetValidationMetrics::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  ++nEvents_;

  edm::Handle<edm::View<reco::Jet>> jetsH;
  edm::Handle<reco::GenJetCollection> genJetsH;

  iEvent.getByToken(jetsToken_, jetsH);
  iEvent.getByToken(genJetsToken_, genJetsH);

  if (nEvents_ == 1) {
    std::cout
      << "jetsTag=" << jetsTag_.encode()
      << " valid=" << jetsH.isValid()
      << " size=" << (jetsH.isValid() ? jetsH->size() : 0) << std::endl;

    std::cout
      << "genJetsTag=" << genJetsTag_.encode()
      << " valid=" << genJetsH.isValid()
      << " size=" << (genJetsH.isValid() ? genJetsH->size() : 0) << std::endl;
      
    std::cout << "[" << moduleDescription().moduleLabel() << "] first-event jets size="
              << (jetsH.isValid() ? jetsH->size() : 0)
              << " genjets size=" << (genJetsH.isValid() ? genJetsH->size() : 0)
              << std::endl;
  
    if (jetsH.isValid() && !jetsH->empty()) {
      const auto& j0 = (*jetsH)[0];
      std::cout << "[" << moduleDescription().moduleLabel() << "] first jet: pt="
                << j0.pt() << " eta=" << j0.eta() << " phi=" << j0.phi()
                << std::endl;
    }
  }

  if (!jetsH.isValid() || !genJetsH.isValid()) return;
  if (jetsH->empty() || genJetsH->empty()) return;

  for (const auto& gjet : *genJetsH) {
    if (gjet.pt() < matchGenPtThreshold_) continue;
    if (std::abs(gjet.eta()) > absEtaMax_) continue;

    ++nGen_;

    int bestIdx = -1;
    double bestDR2 = 1e99;

    for (size_t i = 0; i < jetsH->size(); ++i) {
      const auto& rjet = (*jetsH)[i];
      if (rjet.pt() < recoJetPtThreshold_) continue;

      const double dR2 = reco::deltaR2(gjet.eta(), gjet.phi(), rjet.eta(), rjet.phi());
      if (dR2 < bestDR2) {
        bestDR2 = dR2;
        bestIdx = static_cast<int>(i);
      }
    }

    if (bestIdx < 0) continue;
    if (bestDR2 >= rThreshold_ * rThreshold_) continue;
    if (gjet.pt() == 0.0) continue;

    ++nMatched_;

    const auto& rjet = (*jetsH)[bestIdx];
    const double resp = rjet.pt() / gjet.pt();
    const double dR = std::sqrt(bestDR2);

    sumResp_ += resp;
    sumResp2_ += resp * resp;

    sumDR_ += dR;
    sumDR2_ += dR * dR;

    if (gjet.pt() < 100.0) {
      ++nMatchedLowPt_;
      sumRespLowPt_ += resp;
      sumResp2LowPt_ += resp * resp;
    } else {
      ++nMatchedHighPt_;
      sumRespHighPt_ += resp;
      sumResp2HighPt_ += resp * resp;
    }
  }
}

void JetValidationMetrics::endJob() {
  b_nEvents_ = nEvents_;
  b_nGen_ = nGen_;
  b_nMatched_ = nMatched_;

  if (nMatchedLowPt_ > 0) {
    const double invN = 1.0 / static_cast<double>(nMatchedLowPt_);
    const double mean = sumRespLowPt_ * invN;
    const double var  = (sumResp2LowPt_ * invN) - (mean * mean);
    b_meanRespLowPt_ = mean;
    b_rmsRespLowPt_  = (var > 0.0) ? std::sqrt(var) : 0.0;
  } else {
    b_meanRespLowPt_ = 0.0;
    b_rmsRespLowPt_  = 0.0;
  }

  // High pt (>= 100 GeV)
  if (nMatchedHighPt_ > 0) {
    const double invN = 1.0 / static_cast<double>(nMatchedHighPt_);
    const double mean = sumRespHighPt_ * invN;
    const double var  = (sumResp2HighPt_ * invN) - (mean * mean);
    b_meanRespHighPt_ = mean;
    b_rmsRespHighPt_  = (var > 0.0) ? std::sqrt(var) : 0.0;
  } else {
    b_meanRespHighPt_ = 0.0;
    b_rmsRespHighPt_  = 0.0;
  }

  outTree_->Fill();
}

void JetValidationMetrics::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("jets", edm::InputTag("hltAK4PFPuppiJets"));
  desc.add<edm::InputTag>("genjets", edm::InputTag("ak4GenJetsNoNu"));

  desc.add<double>("recoJetPtThreshold", 30.0);
  desc.add<double>("matchGenPtThreshold", 20.0);
  desc.add<double>("RThreshold", 0.4);
  desc.add<double>("absEtaMax", 6.0);

  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(JetValidationMetrics);
