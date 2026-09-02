#include <algorithm>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "PatternRecognitionbyCLUEstering.h"
#include "TrackstersPCA.h"

using namespace ticl;

template <typename TILES>
PatternRecognitionbyCLUEstering<TILES>::PatternRecognitionbyCLUEstering(const edm::ParameterSet &conf,
                                                                       edm::ConsumesCollector iC)
    : PatternRecognitionAlgoBaseT<TILES>(conf, iC),
      assignmentToken_(iC.consumes<std::vector<int32_t>>(conf.getParameter<edm::InputTag>("tracksterAssignment"))),
      doPidCut_(conf.getParameter<bool>("doPidCut")),
      cutHadProb_(conf.getParameter<double>("cutHadProb")),
      computeLocalTime_(conf.getParameter<bool>("computeLocalTime")),
      usePCACleaning_(conf.getParameter<bool>("usePCACleaning")),
      minNumLayerCluster_(conf.getParameter<int>("minNumLayerCluster")) {}

template <typename TILES>
void PatternRecognitionbyCLUEstering<TILES>::setGeometry(ticlgeom::Tools const &rhtools) {
  this->rhtools_ = &rhtools;
  this->geometryReady_ = true;
}

template <typename TILES>
void PatternRecognitionbyCLUEstering<TILES>::makeTracksters(
    const typename PatternRecognitionAlgoBaseT<TILES>::Inputs &input,
    std::vector<Trackster> &result,
    std::unordered_map<int, std::vector<int>> &seedToTracksterAssociation) {
  if (!this->geometryReady_) {
    throw cms::Exception("PatternRecognitionbyCLUEstering")
        << "setGeometry() must be called before makeTracksters(): geometry is not available.";
  }
  auto const *rhtools = this->rhtools_;

  const auto &assignment = input.ev.get(assignmentToken_);
  if (assignment.size() != input.layerClusters.size()) {
    throw cms::Exception("PatternRecognitionbyCLUEstering")
        << "The trackster assignment has " << assignment.size() << " entries but the merged layer-cluster collection "
        << "has " << input.layerClusters.size()
        << ". The device layer-cluster collections consumed by TrackstersCLUEsteringProducer must be exactly "
        << "the ones merged into the layer-cluster collection used here, in the same order.";
  }

  const auto maxIndex = assignment.empty() ? -1 : *std::max_element(assignment.begin(), assignment.end());
  if (maxIndex < 0) {
    return;
  }
  result.resize(static_cast<size_t>(maxIndex) + 1);

  for (size_t i = 0; i < assignment.size(); ++i) {
    const auto trackster = assignment[i];
    if (trackster < 0 or input.mask[i] == 0.f) {
      continue;
    }
    result[trackster].vertices().push_back(i);
    // CLUE assigns each layer cluster to exactly one trackster, so it is never shared.
    result[trackster].vertex_multiplicity().push_back(1);
  }

  // Same size cut as PatternRecognitionbyCLUE3D (minNumLayerCluster, not inclusive), so that trackster
  // multiplicities are defined identically in the two algorithms.
  result.erase(std::remove_if(result.begin(),
                              result.end(),
                              [this](Trackster const &t) {
                                return static_cast<int>(t.vertices().size()) < minNumLayerCluster_;
                              }),
               result.end());

  // Fills raw energy, barycenter, timing, PCA eigenvalues, eigenvectors and sigmas: without
  // this the tracksters are not usable by the trackster linking and candidate steps.
  const auto limit_em = rhtools->getPositionLayer(rhtools->lastLayerEE(false), false).z();
  ticl::assignPCAtoTracksters(result,
                              input.layerClusters,
                              input.layerClustersTime,
                              limit_em,
                              *rhtools,
                              computeLocalTime_,
                              true,  // energy weighting
                              usePCACleaning_,
                              false);  // never barrel: this plugin is only built for HGCAL tiles

  if (PatternRecognitionAlgoBaseT<TILES>::algo_verbosity_ > VerbosityLevel::Advanced) {
    for (auto const &t : result) {
      edm::LogVerbatim("PatternRecognitionbyCLUEstering")
          << "LCs: " << t.vertices().size() << " Energy: " << t.raw_energy() << " Barycenter: " << t.barycenter();
    }
  }
}

template <typename TILES>
void PatternRecognitionbyCLUEstering<TILES>::filter(
    std::vector<Trackster> &output,
    const std::vector<Trackster> &inTracksters,
    const typename PatternRecognitionAlgoBaseT<TILES>::Inputs &input,
    std::unordered_map<int, std::vector<int>> &seedToTracksterAssociation) {
  auto isHAD = [this](const Trackster &t) -> bool {
    auto const hadProb = t.id_probability(ticl::Trackster::ParticleType::charged_hadron) +
                         t.id_probability(ticl::Trackster::ParticleType::neutral_hadron);
    return hadProb >= cutHadProb_;
  };

  if (doPidCut_) {
    for (auto const &t : inTracksters) {
      if (!isHAD(t)) {
        output.push_back(t);
      }
    }
  } else {
    output = inTracksters;
  }
}

template <typename TILES>
void PatternRecognitionbyCLUEstering<TILES>::fillPSetDescription(edm::ParameterSetDescription &iDesc) {
  iDesc.add<int>("algo_verbosity", 0);
  iDesc.add<edm::InputTag>("tracksterAssignment", edm::InputTag("ticlTrackstersCLUEsteringAssignment"))
      ->setComment("Layer-cluster to trackster assignment produced on device by TrackstersCLUEsteringProducer.");
  iDesc.add<bool>("doPidCut", false);
  iDesc.add<double>("cutHadProb", 0.5);
  iDesc.add<bool>("computeLocalTime", true);
  iDesc.add<bool>("usePCACleaning", true)->setComment("Enable PCA cleaning algorithm");
  iDesc.add<int>("minNumLayerCluster", 2)->setComment("Not inclusive; same definition as CLUE3D");
}

template class ticl::PatternRecognitionbyCLUEstering<TICLLayerTiles>;
