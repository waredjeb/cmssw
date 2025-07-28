#ifndef PhysicsTools_PyTorch_plugins_alpaka_Kernels_h
#define PhysicsTools_PyTorch_plugins_alpaka_Kernels_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/PortableTestObjects/interface/alpaka/TestDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::torchtest {

  /**
   * @brief Utility functions to run simple Alpaka kernels for testing or validation.
   *
   * Simple device-side functionality for modifying and verifying
   * collections of structured SoA data, such as particles, classification outputs, and regressions.
   */
  void fillParticleCollection(Queue &queue, torchportabletest::ParticleCollection &data, float value);
  void assertCombinatorics(Queue &queue, torchportabletest::ParticleCollection &data, float value);
  void assertClassification(Queue &queue, torchportabletest::ClassificationCollection &data);
  void assertRegression(Queue &queue, torchportabletest::RegressionCollection &data);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::torchtest

#endif  // PhysicsTools_PyTorch_plugins_alpaka_Kernels_h
