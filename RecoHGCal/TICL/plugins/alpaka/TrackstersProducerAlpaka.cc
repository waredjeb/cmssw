#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include "CLUEsteringAlgo.h"


namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class TrackstersProducerAlpaka : public global::EDProducer<>{
  public:
    TrackstersProducerAlpaka(edm::ParameterSet const& config) {}
    void produce(edm::StreamID sid, device::Event& event, device::EventSetup const&) const override {
      algo_.printHelloWorld();
      // std::cout << "TRy" << std::endl;
    }

    // implementation of the algorithm
    CLUEsteringAlgo algo_;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(TrackstersProducerAlpaka);
