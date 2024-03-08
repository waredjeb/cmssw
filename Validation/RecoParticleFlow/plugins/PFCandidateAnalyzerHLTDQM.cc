// Producer for particle flow candidates. Plots Eta, Phi, Charge, Pt (log freq, bin)
// for different types of particles described in python/defaults_cfi.py
// It actually uses packedCandidates so that we need only MINIAOD contents to run this DQMAnalyzer.
// note: for pt, log freq is done in this producer, but log freq is done by running
// compare.py
// author: Chosila Sutantawibul, April 23, 2020

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"

#include "TH1F.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>
#include <map>
#include <string>
#include <cstring>

class PFCandidateAnalyzerHLTDQM : public DQMEDAnalyzer {
public:
  explicit PFCandidateAnalyzerHLTDQM(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&) override;

protected:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

private:
  //from config file
  edm::InputTag PFCandTag;
  edm::EDGetTokenT<reco::PFCandidateCollection> PFCandToken;
  std::vector<double> etabins;
  std::map<std::string, MonitorElement*> me;

  std::map<uint32_t, std::string> pdgMap;
};

// constructor
PFCandidateAnalyzerHLTDQM::PFCandidateAnalyzerHLTDQM(const edm::ParameterSet& iConfig) {
  PFCandTag = iConfig.getParameter<edm::InputTag>("PFCandType");
  PFCandToken = consumes<reco::PFCandidateCollection>(PFCandTag);
  etabins = iConfig.getParameter<std::vector<double>>("etabins");

  //create map of pdgId
  std::vector<uint32_t> pdgKeys = iConfig.getParameter<std::vector<uint32_t>>("pdgKeys");
  std::vector<std::string> pdgStrs = iConfig.getParameter<std::vector<std::string>>("pdgStrs");
  for (int i = 0, n = pdgKeys.size(); i < n; i++)
    pdgMap[pdgKeys[i]] = pdgStrs[i];
}

void PFCandidateAnalyzerHLTDQM::bookHistograms(DQMStore::IBooker& booker, edm::Run const&, edm::EventSetup const&) {
  // all candidate
  booker.setCurrentFolder("ParticleFlow/PFCandidate/AllCandidate");

  // for eta binning
  int n = etabins.size() - 1;
  float etabinArray[etabins.size()];
  std::copy(etabins.begin(), etabins.end(), etabinArray);

  //eta has variable bin sizes, use 4th def of TH1F constructor
  TH1F* etaHist = new TH1F("AllCandidateEta", "AllCandidateEta", n, etabinArray);
//  me["AllCandidateEta"] = booker.book1D("AllCandidateEta", etaHist);
  me["AllCandidateEta"] = booker.book1D("AllCandidateEta", "AllCandidateEta", 20, -3.0, 3.0);

  me["AllCandidateLog10Pt"] = booker.book1D("AllCandidateLog10Pt", "AllCandidateLog10Pt", 120, -2, 4);

  //for phi binnings
  double nPhiBins = 73;
  double phiBinWidth = M_PI / (nPhiBins - 1) * 2.;
  me["AllCandidatePhi"] = booker.book1D(
      "AllCandidatePhi", "AllCandidatePhi", nPhiBins, -M_PI - 0.25 * phiBinWidth, +M_PI + 0.75 * phiBinWidth);

  me["AllCandidateCharge"] = booker.book1D("AllCandidateCharge", "AllCandidateCharge", 3, -1.5, 1.5);
  me["AllCandidatePtLow"] = booker.book1D("AllCandidatePtLow", "AllCandidatePtLow", 100, 0., 5.);
  me["AllCandidatePtMid"] = booker.book1D("AllCandidatePtMid", "AllCandidatePtMid", 100, 0., 200.);
  me["AllCandidatePtHigh"] = booker.book1D("AllCandidatePtHigh", "AllCandidatePtHigh", 100, 0., 1000.);
  me["AllCandidateECALEnergyLow"] = booker.book1D("AllCandidateECALEnergyLow", "AllCandidateECALEnergy", 100, 0., 5.);
  me["AllCandidateECALEnergyMid"] = booker.book1D("AllCandidateECALEnergyMid", "AllCandidateECALEnergy", 100, 0., 100.);
  me["AllCandidateECALEnergyHigh"] = booker.book1D("AllCandidateECALEnergyHigh", "AllCandidateECALEnergy", 100, 0., 1000.);
  me["AllCandidateHCALEnergyLow"] = booker.book1D("AllCandidateHCALEnergyLow", "AllCandidateHCALEnergy", 100, 0., 5.);
  me["AllCandidateHCALEnergyMid"] = booker.book1D("AllCandidateHCALEnergyMid", "AllCandidateHCALEnergy", 100, 0., 100.);
  me["AllCandidateHCALEnergyHigh"] = booker.book1D("AllCandidateHCALEnergyHigh", "AllCandidateHCALEnergy", 100, 0., 1000.);
  me["AllCandidateHCALEnergy"] = booker.book1D("AllCandidateHCALEnergy", "AllCandidateHCALEnergy", 100, 0., 1000.);
  me["AllCandidateHOEnergy"] = booker.book1D("AllCandidateHOEnergy", "AllCandidateHOEnergy", 50, 0., 0.1);
  me["AllCandidatePSEnergy"] = booker.book1D("AllCandidatePSEnergy", "AllCandidatePSEnergy", 50, 0., 0.1);
  booker.setCurrentFolder("ParticleFlow/PFCandidate/Undefined");

  me["UndefinedPhi"] = booker.book1D(
      "UndefinedPhi", "UndefinedPhi", nPhiBins, -M_PI - 0.25 * phiBinWidth, +M_PI + 0.75 * phiBinWidth);

  me["UndefinedCharge"] = booker.book1D("UndefinedCharge", "UndefinedCharge", 3, -1.5, 1.5);
  me["UndefinedPtLow"] = booker.book1D("UndefinedPtLow", "UndefinedPtLow", 100, 0., 5.);
  me["UndefinedPtMid"] = booker.book1D("UndefinedPtMid", "UndefinedPtMid", 100, 0., 200.);
  me["UndefinedPtHigh"] = booker.book1D("UndefinedPtHigh", "UndefinedPtHigh", 100, 0., 1000.);
  me["UndefinedECALEnergyLow"] = booker.book1D("UndefinedECALEnergyLow", "UndefinedECALEnergy", 100, 0., 5.);
  me["UndefinedECALEnergyMid"] = booker.book1D("UndefinedECALEnergyMid", "UndefinedECALEnergy", 100, 0., 100.);
  me["UndefinedECALEnergyHigh"] = booker.book1D("UndefinedECALEnergyHigh", "UndefinedECALEnergy", 100, 0., 1000.);
  me["UndefinedHCALEnergyLow"] = booker.book1D("UndefinedHCALEnergyLow", "UndefinedHCALEnergy", 100, 0., 5.);
  me["UndefinedHCALEnergyMid"] = booker.book1D("UndefinedHCALEnergyMid", "UndefinedHCALEnergy", 100, 0., 100.);
  me["UndefinedHCALEnergyHigh"] = booker.book1D("UndefinedHCALEnergyHigh", "UndefinedHCALEnergy", 100, 0., 1000.);
  me["UndefinedHCALEnergy"] = booker.book1D("UndefinedHCALEnergy", "UndefinedHCALEnergy", 100, 0., 1000.);
  me["UndefinedHOEnergy"] = booker.book1D("UndefinedHOEnergy", "UndefinedHOEnergy", 50, 0., 0.1);
  me["UndefinedPSEnergy"] = booker.book1D("UndefinedPSEnergy", "UndefinedPSEnergy", 50, 0., 0.1);

  std::string etaHistName;
  for (auto& pair : pdgMap) {
    booker.setCurrentFolder("ParticleFlow/PFCandidate/" + pair.second);

    //TH1F only takes char*, so have to do conversions for histogram name
    etaHistName = pair.second + "Eta";
    TH1F* etaHist = new TH1F(etaHistName.c_str(), etaHistName.c_str(), n, etabinArray);
//    me[pair.second + "Eta"] = booker.book1D(pair.second + "Eta", etaHist);
    me[pair.second + "Eta"] = booker.book1D(pair.second + "Eta", pair.second + "Eta",  20, -3.0, 3.0);
    me[pair.second + "Log10Pt"] = booker.book1D(pair.second + "Log10Pt", pair.second + "Log10Pt", 120, -2, 4);
    me[pair.second + "Phi"] = booker.book1D(
        pair.second + "Phi", pair.second + "Phi", nPhiBins, -M_PI - 0.25 * phiBinWidth, +M_PI + 0.75 * phiBinWidth);
    me[pair.second + "Charge"] = booker.book1D(pair.second + "Charge", pair.second + "Charge", 3, -1.5, 1.5);
    me[pair.second + "PtLow"] = booker.book1D(pair.second + "PtLow", pair.second + "PtLow", 100, 0., 5.);
    me[pair.second + "PtMid"] = booker.book1D(pair.second + "PtMid", pair.second + "PtMid", 100, 0., 200.);
    me[pair.second + "PtHigh"] = booker.book1D(pair.second + "PtHigh", pair.second + "PtHigh", 100, 0., 1000.);
    me[pair.second + "ECALEnergyLow"] = booker.book1D(pair.second + "ECALEnergyLow", pair.second + "ECALEnergyLow", 100, 0., 5.);
    me[pair.second + "ECALEnergyMid"] = booker.book1D(pair.second + "ECALEnergyMid", pair.second + "ECALEnergyMid", 50, 0., 100.);
    me[pair.second + "ECALEnergyHigh"] = booker.book1D(pair.second + "ECALEnergyHigh", pair.second + "ECALEnergyHigh", 30, 0., 1000.);
    me[pair.second + "HCALEnergyLow"] = booker.book1D(pair.second + "HCALEnergyLow", pair.second + "HCALEnergyLow", 100, 0., 5.);
    me[pair.second + "HCALEnergyMid"] = booker.book1D(pair.second + "HCALEnergyMid", pair.second + "HCALEnergyMid", 50, 0., 100.);
    me[pair.second + "HCALEnergyHigh"] = booker.book1D(pair.second + "HCALEnergyHigh", pair.second + "HCALEnergyHigh", 30, 0., 1000.);
    me[pair.second + "HOEnergy"] = booker.book1D(pair.second + "HOEnergy", pair.second + "HOEnergy", 50, 0., 0.1);
    me[pair.second + "PSEnergy"] = booker.book1D(pair.second + "PSEnergy", pair.second + "PSEnergy", 50, 0., 0.1);
  }
}

void PFCandidateAnalyzerHLTDQM::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //retrieve
  edm::Handle<reco::PFCandidateCollection> pfHandle;
  iEvent.getByToken(PFCandToken, pfHandle);

  if (!pfHandle.isValid()) {
    edm::LogInfo("OutputInfo") << " failed to retrieve data required by ParticleFlow Task";
    edm::LogInfo("OutputInfo") << " ParticleFlow Task cannot continue...!";
    return;
  } else {
    //Analyze
    // Loop Over Particle Flow Candidates

    for (unsigned int i = 0; i < pfHandle->size(); i++) {
      // Fill Histograms for Candidate Methods
      // all candidates
      me["AllCandidateLog10Pt"]->Fill(log10(pfHandle->at(i).pt()));
      me["AllCandidateEta"]->Fill(pfHandle->at(i).eta());
      me["AllCandidatePhi"]->Fill(pfHandle->at(i).phi());
      me["AllCandidateCharge"]->Fill(pfHandle->at(i).charge());
      me["AllCandidatePtLow"]->Fill(pfHandle->at(i).pt());
      me["AllCandidatePtMid"]->Fill(pfHandle->at(i).pt());
      me["AllCandidatePtHigh"]->Fill(pfHandle->at(i).pt());
      me["AllCandidateECALEnergyLow"]->Fill(pfHandle->at(i).rawEcalEnergy());
      me["AllCandidateECALEnergyMid"]->Fill(pfHandle->at(i).rawEcalEnergy());
      me["AllCandidateECALEnergyHigh"]->Fill(pfHandle->at(i).rawEcalEnergy());
      me["AllCandidateHCALEnergyLow"]->Fill(pfHandle->at(i).rawHcalEnergy());
      me["AllCandidateHCALEnergyMid"]->Fill(pfHandle->at(i).rawHcalEnergy());
      me["AllCandidateHCALEnergyHigh"]->Fill(pfHandle->at(i).rawHcalEnergy());
      me["AllCandidateHOEnergy"]->Fill(pfHandle->at(i).rawHoEnergy()); 
      me["AllCandidatePSEnergy"]->Fill(pfHandle->at(i).pS1Energy() + pfHandle->at(i).pS2Energy());

      int pdgId = abs(pfHandle->at(i).pdgId());
      if (pdgMap.find(pdgId) != pdgMap.end()) {
        me[pdgMap[pdgId] + "Log10Pt"]->Fill(log10(pfHandle->at(i).pt()));
        me[pdgMap[pdgId] + "Eta"]->Fill(pfHandle->at(i).eta());
        me[pdgMap[pdgId] + "Phi"]->Fill(pfHandle->at(i).phi());
        me[pdgMap[pdgId] + "Charge"]->Fill(pfHandle->at(i).charge());
        me[pdgMap[pdgId] + "PtLow"]->Fill(pfHandle->at(i).pt());
        me[pdgMap[pdgId] + "PtMid"]->Fill(pfHandle->at(i).pt());
        me[pdgMap[pdgId] + "PtHigh"]->Fill(pfHandle->at(i).pt());
        me[pdgMap[pdgId] + "ECALEnergyLow"]->Fill(pfHandle->at(i).rawEcalEnergy());
        me[pdgMap[pdgId] + "ECALEnergyMid"]->Fill(pfHandle->at(i).rawEcalEnergy());
        me[pdgMap[pdgId] + "ECALEnergyHigh"]->Fill(pfHandle->at(i).rawEcalEnergy());
        me[pdgMap[pdgId] + "HCALEnergyLow"]->Fill(pfHandle->at(i).rawHcalEnergy());
        me[pdgMap[pdgId] + "HCALEnergyMid"]->Fill(pfHandle->at(i).rawHcalEnergy());
        me[pdgMap[pdgId] + "HCALEnergyHigh"]->Fill(pfHandle->at(i).rawHcalEnergy());
        me[pdgMap[pdgId] + "HOEnergy"]->Fill(pfHandle->at(i).rawHoEnergy()); 
        me[pdgMap[pdgId] + "PSEnergy"]->Fill(pfHandle->at(i).pS1Energy() + pfHandle->at(i).pS2Energy());
      }
      else{
        me["UndefinedLog10Pt"]->Fill(log10(pfHandle->at(i).pt()));
        me["UndefinedEta"]->Fill(pfHandle->at(i).eta());
        me["UndefinedPhi"]->Fill(pfHandle->at(i).phi());
        me["UndefinedCharge"]->Fill(pfHandle->at(i).charge());
        me["UndefinedPtLow"]->Fill(pfHandle->at(i).pt());
        me["UndefinedPtMid"]->Fill(pfHandle->at(i).pt());
        me["UndefinedPtHigh"]->Fill(pfHandle->at(i).pt());
        me["UndefinedECALEnergyLow"]->Fill(pfHandle->at(i).rawEcalEnergy());
        me["UndefinedECALEnergyMid"]->Fill(pfHandle->at(i).rawEcalEnergy());
        me["UndefinedECALEnergyHigh"]->Fill(pfHandle->at(i).rawEcalEnergy());
        me["UndefinedHCALEnergyLow"]->Fill(pfHandle->at(i).rawHcalEnergy());
        me["UndefinedHCALEnergyMid"]->Fill(pfHandle->at(i).rawHcalEnergy());
        me["UndefinedHCALEnergyHigh"]->Fill(pfHandle->at(i).rawHcalEnergy());
        me["UndefinedHOEnergy"]->Fill(pfHandle->at(i).rawHoEnergy()); 
        me["UndefinedPSEnergy"]->Fill(pfHandle->at(i).pS1Energy() + pfHandle->at(i).pS2Energy());
      }
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFCandidateAnalyzerHLTDQM);

