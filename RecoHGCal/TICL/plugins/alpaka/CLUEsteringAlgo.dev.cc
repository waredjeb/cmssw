// // Check that ALPAKA_HOST_ONLY is not defined during device compilation:
// #ifdef ALPAKA_HOST_ONLY
// #error ALPAKA_HOST_ONLY defined in device compilation
// #endif

#include <alpaka/alpaka.hpp>
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"
#include <iostream>
#include "CLUEsteringAlgo.h"


namespace ALPAKA_ACCELERATOR_NAMESPACE {
using namespace cms::alpakatools;
void CLUEsteringAlgo::printHelloWorld() const {
    std::cout<<"************Hello World from Alpaka****************"<<std::endl;
}

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE