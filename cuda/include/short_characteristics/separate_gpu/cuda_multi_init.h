#if !defined CUDA_MULTI_INIT_H && defined USE_CUDA
#define CUDA_MULTI_INIT_H

#include "cuda_struct.h"
#include "solvers_struct.h"

namespace cuda {
namespace separate_device {

void InitDirectionsOnMultiDevice(const grid_directions_t &grid_host,
                                 geo::device_host_ptr_t &device_host_ptr,
                                 geo::grid_directions_device_t *&grid_device);

void InitMultiDeviceGrid(int id_dev, const multi_gpu_config_t &gpu_conf, const grid_t &grid_host,
                         const grid_directions_t &grid_dir_host,
                         geo::device_host_ptr_t &device_host_ptr,
                         geo::grid_device_t *&grid_device);

void ClearDirectionsOnMultiDevice(const multi_gpu_config_t &gpu_conf, std::vector<geo::device_host_ptr_t> &device_host_ptrs,
                                  std::vector<geo::grid_directions_device_t *> &grid_devices);

void ClearGridOnMultiDevice(multi_gpu_config_t &gpu_conf, std::vector<geo::device_host_ptr_t> &device_host_ptrs,
                            std::vector<geo::grid_device_t *> &grid_devices);
} // namespace separate_device
} // namespace cuda

#endif //! CUDA_MULTI_INIT_H