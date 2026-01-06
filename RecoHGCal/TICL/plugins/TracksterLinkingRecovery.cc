#include "RecoHGCal/TICL/interface/TracksterLinkingAlgoBase.h"
#include "RecoHGCal/TICL/plugins/TracksterLinkingRecovery.h"

using namespace ticl;

void TracksterLinkingRecovery::linkTracksters(
    const Inputs& input,
    std::vector<Trackster>& resultTracksters,
    std::vector<std::vector<unsigned int>>& linkedResultTracksters,
    std::vector<std::vector<unsigned int>>& linkedTracksterIdToInputTracksterId,
    std::vector<std::vector<float>>& inputTrackstersMasks) {
  const size_t totalSize = input.tracksters.size();
  resultTracksters.reserve(totalSize);

  // Helper to check if a trackster is masked (mask == 0 means masked/skip)
  auto isMasked = [&input, &inputTrackstersMasks](unsigned int globalIdx) {
    const auto& [collIdx, localIdx] = input.tracksters.spanAndLocalIndex(globalIdx);
    return inputTrackstersMasks[collIdx][localIdx] == 0.f;
  };

  // Pass through all non-masked tracksters
  for (size_t globalIdx = 0; globalIdx < totalSize; ++globalIdx) {
    // Skip masked tracksters
    if (isMasked(globalIdx))
      continue;

    resultTracksters.push_back(input.tracksters[globalIdx]);
    linkedResultTracksters.emplace_back(std::initializer_list<unsigned int>{
        static_cast<unsigned int>(resultTracksters.size() - 1)});
    linkedTracksterIdToInputTracksterId.emplace_back(std::initializer_list<unsigned int>{
        static_cast<unsigned int>(globalIdx)});
    // Mark this trackster as used in the input mask
    const auto& [collectionIdx, localIdx] = input.tracksters.spanAndLocalIndex(globalIdx);
    inputTrackstersMasks[collectionIdx][localIdx] = 0.f;
  }
}
