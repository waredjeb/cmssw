#ifndef PhysicsTools_PyTorch_interface_Converter_h
#define PhysicsTools_PyTorch_interface_Converter_h

#include <torch/torch.h>

#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "PhysicsTools/PyTorch/interface/SoAMetadata.h"

namespace cms::torch::alpaka {

  // Metadata to run model with input SOA and fill output SOA.
  template <typename SOA_Input, typename SOA_Output>
  class ModelMetadata {
  public:
    SoAMetadata<SOA_Input> input;
    SoAMetadata<SOA_Output> output;

    // Used in AOT model class to correctly choose multi or single output conversion
    // Default value true, as single value can be parsed with multi output
    bool multi_output;

    ModelMetadata(const SoAMetadata<SOA_Input>& input_,
                  const SoAMetadata<SOA_Output>& output_,
                  bool multi_output_ = true)
        : input(input_), output(output_), multi_output(multi_output_) {}
  };

  // Static class to wrap raw SOA pointer in tensor object without copying.
  class Converter {
  public:
    // Calculate size and stride of data store based on InputMetadata and return list of IValue, which is parent class of torch::tensor.
    template <typename SOA_Input, typename SOA_Output>
    static std::vector<::torch::IValue> convert_input(const ModelMetadata<SOA_Input, SOA_Output>& metadata,
                                                      ::torch::Device device) {
      std::vector<::torch::IValue> tensors(metadata.input.nBlocks);
      for (int i = 0; i < metadata.input.nBlocks; i++) {
        assert(reinterpret_cast<intptr_t>(metadata.input[metadata.input.order[i]].ptr_) % SOA_Input::alignment == 0);
        tensors[i] = Converter::array_to_tensor(device, metadata.input[metadata.input.order[i]]);
      }
      return tensors;
    }

    // AOT specific implementation, as model expects vector of torch::Tensor not torch::IValue
    template <typename SOA_Input, typename SOA_Output>
    static std::vector<::torch::Tensor> convert_input_tensor(const ModelMetadata<SOA_Input, SOA_Output>& metadata,
                                                             ::torch::Device device) {
      std::vector<::torch::Tensor> tensors(metadata.input.nBlocks);
      for (int i = 0; i < metadata.input.nBlocks; i++) {
        assert(reinterpret_cast<intptr_t>(metadata.input[metadata.input.order[i]].ptr_) % SOA_Input::alignment == 0);
        tensors[i] = Converter::array_to_tensor(device, metadata.input[metadata.input.order[i]]);
      }
      return tensors;
    }

    // Calculate size and stride of data store based on OutputMetadata and return single output tensor
    template <typename SOA_Input, typename SOA_Output>
    static ::torch::Tensor convert_output(const ModelMetadata<SOA_Input, SOA_Output>& metadata,
                                          ::torch::Device device) {
      assert(reinterpret_cast<intptr_t>(metadata.output[metadata.output.order[0]].ptr_) %
                 metadata.output[metadata.output.order[0]].alignment_ ==
             0);
      return Converter::array_to_tensor(device, metadata.output[metadata.output.order[0]]);
    }

    // Calculate size and stride of data store based on OutputMetadata and fill SoA with tensor values
    template <typename SOA_Input, typename SOA_Output>
    static void convert_output(const std::vector<::torch::IValue>& tensors,
                               const ModelMetadata<SOA_Input, SOA_Output>& metadata,
                               ::torch::Device device) {
      for (int i = 0; i < metadata.output.nBlocks; i++) {
        // Only tensors are currenlty supported for conversion
        if (tensors[i].isTensor()) {
          assert(reinterpret_cast<intptr_t>(metadata.output[metadata.output.order[i]].ptr_) % SOA_Output::alignment ==
                 0);
          Converter::array_to_tensor(device, metadata.output[metadata.output.order[i]]) = tensors[i].toTensor();
        }
      }
    }

    // AOT specific implementation, as return type is torch::Tensor not torch::IValue
    template <typename SOA_Input, typename SOA_Output>
    static void convert_output(const std::vector<::torch::Tensor>& tensors,
                               const ModelMetadata<SOA_Input, SOA_Output>& metadata,
                               ::torch::Device device) {
      for (int i = 0; i < metadata.output.nBlocks; i++) {
        assert(reinterpret_cast<intptr_t>(metadata.output[metadata.output.order[i]].ptr_) % SOA_Output::alignment == 0);
        Converter::array_to_tensor(device, metadata.output[metadata.output.order[i]]) = tensors[i];
      }
    }

  private:
    // Wrap raw pointer by torch::Tensor based on type, size and stride.
    static ::torch::Tensor array_to_tensor(::torch::Device device, const Block& block) {
      auto options = ::torch::TensorOptions()
                         .dtype(block.type_)
#ifdef ALPAKA_ACC_GPU_CUDA_ENABLED
                         .device(device)
#endif
                         .pinned_memory(true);
      return ::torch::from_blob(block.ptr_, block.size_, block.stride_, options);
    }
  };

}  // namespace cms::torch::alpaka

#endif  // PhysicsTools_PyTorch_interface_Converter_h
