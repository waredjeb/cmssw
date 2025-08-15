#include <alpaka/alpaka.hpp>
#include <torch/script.h>
#include <torch/torch.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#ifdef ALPAKA_ACC_GPU_CUDA_ENABLED
#include <c10/cuda/CUDAStream.h>
#include <cuda_runtime.h>
#endif

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "DataFormats/Portable/interface/PortableCollection.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "PhysicsTools/PyTorch/interface/AlpakaConfig.h"
#include "PhysicsTools/PyTorch/interface/Converter.h"
#include "PhysicsTools/PyTorch/interface/Model.h"

#include "PhysicsTools/PyTorch/test/testBase.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using JitModel = cms::torch::alpaka::Model<cms::torch::alpaka::CompilationType::kJustInTime>;

  // Input SOA
  GENERATE_SOA_LAYOUT(SoAPositionTemplate, SOA_COLUMN(float, x), SOA_COLUMN(float, y), SOA_COLUMN(float, z))
  GENERATE_SOA_LAYOUT(SoARotationTemplate, SOA_COLUMN(float, eta), SOA_COLUMN(float, phi))
  GENERATE_SOA_LAYOUT(SoAClusterTemplate, SOA_EIGEN_COLUMN(Eigen::Vector3d, position), SOA_SCALAR(int, count))

  using SoAPosition = SoAPositionTemplate<>;
  using SoARotation = SoARotationTemplate<>;
  using SoACluster = SoAClusterTemplate<>;

  // Output SOA
  GENERATE_SOA_LAYOUT(SoAResultTemplate, SOA_COLUMN(float, x), SOA_COLUMN(float, y))
  GENERATE_SOA_LAYOUT(SoAResultClusterTemplate, SOA_EIGEN_COLUMN(Eigen::Matrix3f, outer))

  using SoAResult = SoAResultTemplate<>;
  using SoAResultCluster = SoAResultClusterTemplate<>;

  class testMultiSoA : public testBasePyTorch {
    CPPUNIT_TEST_SUITE(testMultiSoA);
    CPPUNIT_TEST(test);
    CPPUNIT_TEST_SUITE_END();

  public:
    std::string pyScript() const override;
    void test() override;
  };

  std::string testMultiSoA::pyScript() const { return "create_multi_input_output_model.py"; }

  CPPUNIT_TEST_SUITE_REGISTRATION(testMultiSoA);

  class FillKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<::alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  PortableCollection<SoAPosition, Device>::View positionView,
                                  PortableCollection<SoARotation, Device>::View rotationView,
                                  PortableCollection<SoACluster, Device>::View clusterView) const {
      if (cms::alpakatools::once_per_grid(acc)) {
        clusterView.count() = clusterView.metadata().size();
      }

      for (int32_t i : cms::alpakatools::uniform_elements(acc, positionView.metadata().size())) {
        positionView.x()[i] = i;
        positionView.y()[i] = i + 1;
        positionView.z()[i] = i + 2;

        rotationView.eta()[i] = i;
        rotationView.phi()[i] = i * 2;
      }

      for (int32_t i : cms::alpakatools::uniform_elements(acc, clusterView.metadata().size())) {
        clusterView[i].position()(0) = i;
        clusterView[i].position()(1) = i + 1;
        clusterView[i].position()(2) = i + 2;
      }
    }
  };

  class TestVerifyKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc, PortableCollection<SoAResult, Device>::View view) const {
      float result_check[4][2] = {{6, 6}, {9, 12}, {18, 24}, {33, 42}};
      for (uint32_t i : cms::alpakatools::uniform_elements(acc, view.metadata().size())) {
        ALPAKA_ASSERT_ACC(view.x()[i] - result_check[i][0] < 1.0e-05);
        ALPAKA_ASSERT_ACC(view.x()[i] - result_check[i][0] > -1.0e-05);
        ALPAKA_ASSERT_ACC(view.y()[i] - result_check[i][1] < 1.0e-05);
        ALPAKA_ASSERT_ACC(view.y()[i] - result_check[i][1] > -1.0e-05);
      }
    }
  };

  void fill(Queue& queue,
            PortableCollection<SoAPosition, Device>& positionCollection,
            PortableCollection<SoARotation, Device>& rotationCollection,
            PortableCollection<SoACluster, Device>& clusterCollection) {
    uint32_t items = 64;
    uint32_t groups = cms::alpakatools::divide_up_by(12, items);
    auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(groups, items);
    alpaka::exec<Acc1D>(
        queue, workDiv, FillKernel{}, positionCollection.view(), rotationCollection.view(), clusterCollection.view());
  }

  void check(Queue& queue, PortableCollection<SoAResult, Device>& collection) {
    uint32_t items = 64;
    uint32_t groups = cms::alpakatools::divide_up_by(collection->metadata().size(), items);
    auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(groups, items);
    alpaka::exec<Acc1D>(queue, workDiv, TestVerifyKernel{}, collection.view());
  }

  void testMultiSoA::test() {
    Platform platform;
    std::vector<Device> alpakaDevices = ::alpaka::getDevs(platform);
    const auto& alpakaHost = ::alpaka::getDevByIdx(alpaka_common::PlatformHost(), 0u);
    CPPUNIT_ASSERT(alpakaDevices.size());
    const auto& alpakaDevice = alpakaDevices[0];
    Queue queue{alpakaDevice};
    torch::Device torchDevice(cms::torch::alpaka::kTorchDeviceType);

    // Number of elements
    const std::size_t element_count = 4;
    const std::size_t cluster_size = 6;

    // Create and fill needed portable collections
    PortableCollection<SoAPosition, Device> positionCollection(element_count, alpakaDevice);
    PortableCollection<SoARotation, Device> rotationCollection(element_count, alpakaDevice);
    PortableCollection<SoACluster, Device> clusterCollection(cluster_size, alpakaDevice);

    PortableCollection<SoAResult, Device> resultCollection(element_count, alpakaDevice);
    // PortableCollection<SoAResultCluster, Device> resultClusterCollection(cluster_size, alpakaDevice);
    fill(queue, positionCollection, rotationCollection, clusterCollection);
    alpaka::wait(queue);

    std::string model_path = dataPath_ + "/multi_input_output_model.pt";
    auto model = JitModel(model_path);
    model.to(queue);

    // Create SoA Metadata
    cms::torch::alpaka::SoAMetadata<SoAPosition> input(element_count);
    auto positionRecords = positionCollection.view().records();
    auto rotationRecords = rotationCollection.view().records();
    auto clusterRecords = clusterCollection.view().records();
    input.append_block("position", positionRecords.x(), positionRecords.y(), positionRecords.z());
    input.append_block<SoARotation>("rotation", element_count, rotationRecords.eta(), rotationRecords.phi());
    input.append_block<SoACluster>("cluster", cluster_size, clusterRecords.position());
    input.append_block<SoACluster>("clusterCount", element_count, clusterRecords.count());

    cms::torch::alpaka::SoAMetadata<SoAResult> output(element_count);
    auto resultView = resultCollection.view().records();
    // auto resultClusterView = resultClusterCollection.view().records();
    output.append_block("result", resultView.x(), resultView.y());

    // TODO: Make multi output
    // output.append_block<SoAResultCluster>("cluster", cluster_size, resultClusterView.outer());

    cms::torch::alpaka::ModelMetadata metadata(input, output);

    // Call function to build tensor and run model
    model.forward(metadata);
    alpaka::wait(queue);

    check(queue, resultCollection);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE