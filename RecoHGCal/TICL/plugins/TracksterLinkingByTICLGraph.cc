#include "RecoHGCal/TICL/interface/TracksterLinkingAlgoBase.h"
#include "fastjet/ClusterSequence.hh"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoHGCal/TICL/plugins/TracksterLinkingByTICLGraph.h"

using namespace ticl;

void TracksterLinkingByTICLGraph::linkTracksters(
    const Inputs& input,
    std::vector<Trackster>& resultTracksters,
    std::vector<std::vector<unsigned int>>& linkedResultTracksters,
    std::vector<std::vector<unsigned int>>& linkedTracksterIdToInputTracksterId) {

  auto& ticlGraph = input.ticlGraph;

  auto const subComponents = ticlGraph.findSubComponents();
  linkedTracksterIdToInputTracksterId.resize(subComponents.size());
  // Link tracksters based on which ones are components of the same jet
  for (unsigned int i = 0; i < subComponents.size(); ++i) {
    const auto& subC = subComponents[i];

    std::vector<unsigned int> linkedTracksters;
    Trackster outTrackster;
    if (!subC.empty()) {
      // Check if a trackster is a component of the current jet
      for (const auto& tracksterIndex : subC) {
        linkedTracksterIdToInputTracksterId[i].push_back(tracksterIndex);
        outTrackster.mergeTracksters(input.tracksters[tracksterIndex]);
      }
      linkedTracksters.push_back(resultTracksters.size());
      resultTracksters.push_back(outTrackster);
      // Store the linked tracksters
      linkedResultTracksters.push_back(linkedTracksters);
    }
  }
}
