// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include "PhysicsTools/PyTorch/plugins/alpaka/Kernels.h"

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::torchtest {

  using namespace cms::alpakatools;

  /**
   * @brief Fill all values in a particle collection with a specified constant.
   * 
   * For debugging and unit testing.
   *
   * @param queue Alpaka execution queue.
   * @param data Particle collection to be modified.
   * @param value Constant value to fill the collection with.
   */
  void fillParticleCollection(Queue &queue, torchportabletest::ParticleCollection &data, float value) {
    uint32_t threads_per_block = 512;
    uint32_t blocks_per_grid = divide_up_by(data.view().metadata().size(), threads_per_block);
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);
    alpaka::exec<Acc1D>(
        queue,
        grid,
        [] ALPAKA_FN_ACC(Acc1D const &acc, torchportabletest::ParticleCollection::View data, float value) {
          for (auto tid : uniform_elements(acc, data.metadata().size())) {
            data.pt()[tid] = value;
            data.phi()[tid] = value;
            data.eta()[tid] = value;
          }
        },
        data.view(),
        value);
  }

  /**
   * @brief Assert that the particle collection obeys certain combinatoric relationships.
   *
   * Used in test scenarios to verify data layout or transformation logic.
   *
   * @param queue Alpaka execution queue.
   * @param data Particle collection to check.
   * @param value Reference value for validation.
   */
  void assertCombinatorics(Queue &queue, torchportabletest::ParticleCollection &data, float value) {
    uint32_t threads_per_block = 512;
    uint32_t blocks_per_grid = divide_up_by(data.view().metadata().size(), threads_per_block);
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);
    alpaka::exec<Acc1D>(
        queue,
        grid,
        [] ALPAKA_FN_ACC(Acc1D const &acc, torchportabletest::ParticleCollection::View data, float value) {
          for (auto tid : uniform_elements(acc, data.metadata().size())) {
            ALPAKA_ASSERT_ACC(data.pt()[tid] == value);
            ALPAKA_ASSERT_ACC(data.phi()[tid] == value);
            ALPAKA_ASSERT_ACC(data.eta()[tid] == value);
          }
        },
        data.view(),
        value);
  }

  /**
   * @brief Validate classification model outputs.
   *
   * Checks whether the classification outputs match expected format or values.
   * Used for debugging and integration testing of model inference.
   *
   * @param queue Alpaka execution queue.
   * @param data Classification output collection.
   */
  void assertClassification(Queue &queue, torchportabletest::ClassificationCollection &data) {
    uint32_t threads_per_block = 512;
    uint32_t blocks_per_grid = divide_up_by(data.view().metadata().size(), threads_per_block);
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);
    alpaka::exec<Acc1D>(
        queue,
        grid,
        [] ALPAKA_FN_ACC(Acc1D const &acc, torchportabletest::ClassificationCollection::View data) {
          for (auto tid : uniform_elements(acc, data.metadata().size())) {
            ALPAKA_ASSERT_ACC(data.c1()[tid] == 0.5f);
            ALPAKA_ASSERT_ACC(data.c2()[tid] == 0.5f);
          }
        },
        data.view());
  }

  /**
   * @brief Validate regression model outputs.
   *
   * Similar to classification checks, this is used for asserting output correctness in regression tasks.
   *
   * @param queue Alpaka execution queue.
   * @param data Regression output collection.
   */
  void assertRegression(Queue &queue, torchportabletest::RegressionCollection &data) {
    uint32_t threads_per_block = 512;
    uint32_t blocks_per_grid = divide_up_by(data.view().metadata().size(), threads_per_block);
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);
    alpaka::exec<Acc1D>(
        queue,
        grid,
        [] ALPAKA_FN_ACC(Acc1D const &acc, torchportabletest::RegressionCollection::View data) {
          for (auto tid : uniform_elements(acc, data.metadata().size())) {
            ALPAKA_ASSERT_ACC(data.reco_pt()[tid] == 0.5f);
          }
        },
        data.view());
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::torchtest
