#ifndef RecoHGCal_TICL_TracksterCleaningByBeta_H
#define RecoHGCal_TICL_TracksterCleaningByBeta_H

#include <string>
#include "RecoHGCal/TICL/interface/TracksterCleaningAlgoBase.h"

namespace ticl {

class TracksterCleaningByBeta final : public TracksterCleaningAlgoBase {
public:
  TracksterCleaningByBeta(const edm::ParameterSet& conf, edm::ConsumesCollector);

  void cleanTracksters(const Inputs& input,
                       std::vector<ticl::Trackster>& outTracksters,
                       std::vector<std::vector<unsigned int>>& outMap) const override;

  void initialize() override {}
  void setEvent(edm::Event&, edm::EventSetup const&) override {}

  static void fillPSetDescription(edm::ParameterSetDescription& desc) {
    TracksterCleaningAlgoBase::fillPSetDescription(desc);
    desc.add<double>("betaContamMin", 1.12);
    desc.add<double>("R0", 0.10);
    desc.add<bool>("useRawEnergy", true);
    desc.add<double>("epsE", 1e-6);
    desc.add<double>("epsDR", 1e-6);
    desc.add<bool>("weightMode", true);
    desc.add<bool>("emitDroppedAsStandalone", false);
    desc.add<double>("zAbsCut", 25.0);
    desc.add<double>("tAbsCut", 0.15);
    desc.add<double>("sigmaZ", 12.5);
    desc.add<double>("sigmaT", 0.08);
    desc.add<double>("sigmaDR", 0.08);
    desc.add<double>("zPower", 1.5);
    desc.add<double>("tPower", 0.5);
    desc.add<double>("drPower", 0.5);
    desc.add<double>("wmin", 1e-3);

    desc.add<bool>("doPruning", true);      
    desc.add<double>("pruneWmin", 1e-2);              
    desc.add<bool>("pruneUseSeparateKernels", false);  
    desc.add<double>("sigmaZ_prune", 12.5);
    desc.add<double>("sigmaT_prune", 0.08);
    desc.add<double>("sigmaDR_prune", 0.08);
    desc.add<double>("zPower_prune", 1.5);
    desc.add<double>("tPower_prune", 0.5);
    desc.add<double>("drPower_prune", 0.5);
  }

private:
  double betaContamMin_, R0_, epsE_, epsDR_;
  bool   useRawEnergy_, emitDroppedAsStandalone_;
  bool   weightMode_;
  double zAbsCut_, tAbsCut_;
  double sigmaZ_, sigmaT_, sigmaDR_;
  double zPower_, tPower_, drPower_;
  double wmin_;

  // pruning
  bool   doPruning_; 
  double pruneWmin_;   
  bool   pruneUseSeparateKernels_;
  double sigmaZ_prune_, sigmaT_prune_, sigmaDR_prune_;
  double zPower_prune_, tPower_prune_, drPower_prune_;


  inline double linkEnergy_(const ticl::Trackster& trk) const {
    return useRawEnergy_ ? trk.raw_energy() : trk.regressed_energy();
  }
  inline void setLinkRawEnergy_(ticl::Trackster& trk, double e) const {
    trk.setRawEnergy(static_cast<float>(e));
  }
};

}  // namespace ticl

#endif
