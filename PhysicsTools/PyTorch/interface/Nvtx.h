#ifndef PhysicsTools_PyTorch_interface_Nvtx_h
#define PhysicsTools_PyTorch_interface_Nvtx_h

#ifdef ALPAKA_ACC_GPU_CUDA_ENABLED
#include <nvtx3/nvToolsExt.h>
#endif

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  /**
 * @class NvtxScopedRange
 * @brief Helper class for NVTX profiling.
 *
 * Exposes a simple interface to create and manage NVTX ranges.
 * Automatically ends the range when the object goes out of scope.
 *
 * Only enabled on CUDA backend via NVTX.
 */
  class NvtxScopedRange {
  public:
    explicit NvtxScopedRange(const char* msg) {
#ifdef ALPAKA_ACC_GPU_CUDA_ENABLED
      id_ = nvtxRangeStartA(msg);
#endif
    }

    void end() {
#ifdef ALPAKA_ACC_GPU_CUDA_ENABLED
      if (active_) {
        active_ = false;
        nvtxRangeEnd(id_);
      }
#endif
    }

    ~NvtxScopedRange() { end(); }

  private:
#ifdef ALPAKA_ACC_GPU_CUDA_ENABLED
    nvtxRangeId_t id_;
    bool active_ = true;
#endif
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // PhysicsTools_PyTorch_interface_Nvtx_h
