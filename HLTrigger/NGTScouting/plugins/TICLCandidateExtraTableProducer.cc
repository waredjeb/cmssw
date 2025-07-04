#include <algorithm>
#include <numeric>
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TICLCandidateExtraTableProducer : public edm::global::EDProducer<> {
public:
  TICLCandidateExtraTableProducer(const edm::ParameterSet& cfg)
      : tableName_(cfg.getParameter<std::string>("tableName")),
        skipNonExistingSrc_(cfg.getParameter<bool>("skipNonExistingSrc")),
        ticlCandidate_token_(mayConsume<std::vector<TICLCandidate>>(cfg.getParameter<edm::InputTag>("candidates"))),
        precision_(cfg.getParameter<int>("precision")) {
    produces<nanoaod::FlatTable>(tableName_);         // main table
    produces<nanoaod::FlatTable>("tracksterIndices"); // child table
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("tableName", "hltTiclCandidatesTable")
        ->setComment("Name of the main TICLCandidate table");
    desc.add<bool>("skipNonExistingSrc", false)
        ->setComment("Skip if missing input");
    desc.add<edm::InputTag>("candidates", edm::InputTag("hltTiclCandidates"));
    desc.add<int>("precision", 7);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  void produce(edm::StreamID id, edm::Event& event, const edm::EventSetup& setup) const override {
    const auto& ticlCandidatesHandle = event.getHandle(ticlCandidate_token_);
    if (!ticlCandidatesHandle.isValid()) {
      if (skipNonExistingSrc_) return;
      throw cms::Exception("MissingProduct") << "TICLCandidates not found!\n";
    }

    const auto& ticlCandidates = *ticlCandidatesHandle;
    const size_t nCandidates = ticlCandidates.size();

    // One flat vector for ALL trackster indices
    std::vector<int> flatTracksterIndices;
    // One count per TICLCandidate
    std::vector<uint32_t> nTrackstersPerCandidate;
    flatTracksterIndices.reserve(nCandidates * 5);  // guess: avg 5 per candidate
    nTrackstersPerCandidate.reserve(nCandidates);

    for (const auto& candidate : ticlCandidates) {
      const auto& tracksters = candidate.tracksters();
      uint32_t count = 0;
      for (const auto& trackster : tracksters) {
        if (trackster.isNonnull()) {
          flatTracksterIndices.push_back(trackster.key());
          ++count;
        }
      }
      nTrackstersPerCandidate.push_back(count);
    }

    // === Main table ===
    auto mainTable = std::make_unique<nanoaod::FlatTable>(nCandidates, tableName_, /*singleton=*/false, /*extension=*/true);
    mainTable->setDoc("TICLCandidates with associated Tracksters");
    mainTable->addColumn<uint32_t>(
      "nTracksters",
      nTrackstersPerCandidate,
      "Number of Tracksters associated to this TICLCandidate"
    );

    // === Child table ===
    auto childTable = std::make_unique<nanoaod::FlatTable>(
      flatTracksterIndices.size(), "tracksterIndices", /*singleton=*/false, /*extension=*/true);
    childTable->setDoc("Flattened indices of all Tracksters linked to TICLCandidates");
    childTable->addColumn<int>(
      "tracksterIndex",
      flatTracksterIndices,
      "Index in the Trackster collection for each Trackster in TICLCandidates",
      precision_
    );

    // Put both tables
    event.put(std::move(mainTable), tableName_);
    event.put(std::move(childTable), "tracksterIndices");
  }

private:
  const std::string tableName_;
  const bool skipNonExistingSrc_;
  const edm::EDGetTokenT<std::vector<TICLCandidate>> ticlCandidate_token_;
  const unsigned int precision_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TICLCandidateExtraTableProducer);

