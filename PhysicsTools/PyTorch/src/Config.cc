#include "PhysicsTools/PyTorch/interface/Config.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace cms::torch {

  ::torch::jit::script::Module load(const std::string &model_path) {
    ::torch::jit::script::Module model;
    try {
      model = ::torch::jit::load(model_path);
    } catch (const c10::Error &e) {
      cms::Exception ex("ModelLoadingError");
      ex << "Error loading the model from path: " << model_path << "\n"
         << "Details: " << e.what();
      ex.addContext("Calling cms::torch::load(const std::string&)");
      throw ex;
    }
    return model;
  }

}  // namespace cms::torch