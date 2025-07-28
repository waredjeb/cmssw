#ifndef PhysicsTools_PyTorch_interface_Config_h
#define PhysicsTools_PyTorch_interface_Config_h

#ifdef ClassDef
#undef ClassDef
#endif

#include <torch/script.h>
#include <torch/torch.h>

namespace cms::torch {

  /**
   * @brief Loads a TorchScript model.
   *
   * This function wraps `torch::jit::load` to load a TorchScript model from a specified path.
   * In case of failure, it throws a CMS-specific `cms::Exception` with detailed context and error information.
   *
   * @param model_path The file path to the TorchScript model (.pt file).
   * @return A loaded `torch::jit::script::Module` ready for inference or further manipulation.
   *
   * @throws cms::Exception If loading fails due to file issues, format errors, or internal TorchScript problems.
   *
   * @note This function is intended for model loading in CMSSW environments, providing
   *       integration with the framework's exception handling and logging facilities.
   */
  ::torch::jit::script::Module load(const std::string &model_path);

}  // namespace cms::torch

#endif  // PhysicsTools_PyTorch_interface_Config_h
