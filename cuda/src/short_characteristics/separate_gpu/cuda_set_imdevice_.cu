#ifdef USE_CUDA
#include "cuda_illum_sep_param.h"

namespace cuda_sep = cuda::separate_device::kernel;
__global__ void cuda_sep::SetImDevice(geo::grid_device_t *grid,
                                      const IdType loc_size_gpu,
                                      const IdType shift_gpu,
                                      const IdType loc_size_params,
                                      const IdType shift_params) {

  grid->loc_size_gpu = loc_size_gpu;
  grid->shift_gpu = shift_gpu;
  grid->loc_size_params = loc_size_params;
  grid->shift_params = shift_params;
}
#endif