#include "cuda_init_mem.h"
#include "cuda_memory.h"
#include "cuda_multi_init.h"
#include "cuda_multi_interface.h"

#include "mpi_shifts.h"

#ifdef SEPARATE_GPU
cuda::multi_gpu_config_t gpu_config;

std::vector<cuda::geo::grid_directions_device_t *> grid_dir_deviceN;
std::vector<cuda::geo::grid_device_t *> grid_deviceN;
std::vector<cuda::geo::device_host_ptr_t> device_host_ptrN;

static void init_params_config(cuda::multi_gpu_config_t &gpu_conf, const grid_t &grid_host) {

  gpu_conf.size_params.resize(gpu_conf.GPU_N, 0);
  gpu_conf.disp_params.resize(gpu_conf.GPU_N, 0);

  const IdType cpu_size = grid_host.loc_size;
  const IdType cpu_shift = grid_host.loc_shift;

  for (size_t id_dev = 0; id_dev < gpu_conf.GPU_N; id_dev++) {

    const IdType gpu_size = gpu_conf.size[id_dev];
    const IdType gpu_disp = gpu_conf.disp[id_dev];

    //левая граница на узле больше правой на карте
    if (cpu_shift > gpu_disp + gpu_size) {
      continue;
    }

    //правая граница на узле меньше левой на карте
    if (cpu_shift + cpu_size < gpu_disp) {
      continue;
    }

    //здесь уже точно есть пересечение
    if (cpu_shift < gpu_disp) {
      int left = 0;
      int right = std::min(cpu_shift + cpu_size, gpu_disp + gpu_size) - gpu_disp;

      gpu_conf.disp_params[id_dev] = left;
      gpu_conf.size_params[id_dev] = right - left;
    } else {
      int left = cpu_shift - gpu_disp;
      int right =
          std::min(std::min(cpu_shift + cpu_size, gpu_disp + gpu_size), gpu_size);

      gpu_conf.disp_params[id_dev] = left;
      gpu_conf.size_params[id_dev] = right - left;
    }
  }
#ifdef LOG_SEP_CUDA
  for (size_t id_dev = 0; id_dev < gpu_conf.GPU_N; id_dev++) {

    const IdType gpu_size = gpu_conf.size[id_dev];
    const IdType gpu_disp = gpu_conf.disp[id_dev];
    log_sep_cuda("[%ld, %ld]gpu_cfg %ld %ld gpu_params %ld %ld\n", cpu_shift,
                 cpu_size, gpu_disp, gpu_size, gpu_conf.disp_params[id_dev], gpu_conf.size_params[id_dev]);
  }
#endif
}

int cuda::interface::separate_device::InitDevice(const grid_directions_t &grid_dir_host, grid_t &grid_host) {

  CUDA_CALL_FUNC(cudaGetDeviceCount, &gpu_config.GPU_N);
#ifdef SINGLE_GPU
  gpu_config.GPU_N = GPU_DIV_PARAM;
#endif
  GetSend(gpu_config.GPU_N, grid_host.size, gpu_config.size);
  GetDisp(gpu_config.GPU_N, grid_host.size, gpu_config.disp);
  init_params_config(gpu_config, grid_host);

  int _dev_n = gpu_config.GPU_N; //реальное число карт
#ifdef SINGLE_GPU
  _dev_n = 1;
#endif

  device_host_ptrN.resize(_dev_n);
  grid_dir_deviceN.resize(_dev_n);
  grid_deviceN.resize(_dev_n);

  for (int dev_id = 0; dev_id < _dev_n; dev_id++) {
    CUDA_CALL_FUNC(cudaSetDevice, dev_id);
    cuda::separate_device::InitDirectionsOnMultiDevice(grid_dir_host, device_host_ptrN[dev_id], grid_dir_deviceN[dev_id]);
    cuda::separate_device::InitMultiDeviceGrid(dev_id, gpu_config, grid_host, grid_dir_host, device_host_ptrN[dev_id], grid_deviceN[dev_id]);
  }

  /// \todo это отдельно, т.к. относится к инициализации хоста

  IdType frq_cnt = 1;
#ifdef SAVE_FULL_SPECTRUM
  frq_cnt = grid_host.size_frq;
#endif

  mem_protected::MallocHost((grid_dir_host.size * grid_host.size * frq_cnt * sizeof(Type)), &grid_host.Illum);
  mem_protected::MallocHost(((grid_dir_host.loc_size) * grid_host.size * frq_cnt * sizeof(Type)), &grid_host.scattering);

  WRITE_LOG("grid_host.loc_size =%lu %lu %lu %lu\n", grid_host.loc_size, grid_host.loc_shift, grid_host.size, grid_host.size_face);

#ifdef ON_FULL_ILLUM_ARRAYS
  if (grid_dir_host.loc_shift == 0 || (grid_host.size != grid_host.loc_size)) // или нулевой узел, или данные разделены
  {
    mem_protected::MallocHost((grid_host.loc_size * sizeof(Type)), &grid_host.energy);
    mem_protected::MallocHost((grid_host.loc_size * sizeof(Vector3)), &grid_host.stream);
    mem_protected::MallocHost((grid_host.loc_size * sizeof(Matrix3)), &grid_host.impuls);
  }
#endif

  return e_completion_success;
}

void cuda::interface::separate_device::ClearDevice() {

  cuda::separate_device::ClearDirectionsOnMultiDevice(gpu_config, device_host_ptrN, grid_dir_deviceN);
  cuda::separate_device::ClearGridOnMultiDevice(gpu_config, device_host_ptrN, grid_deviceN);

  // CUDA_CALL_FUNC(cudaDeviceReset);

  WRITE_LOG("Free device arrays\n");
}
void cuda::interface::separate_device::ClearHost(grid_t &grid_host) {

  mem_protected::FreeMemHost(grid_host.Illum);
  mem_protected::FreeMemHost(grid_host.scattering);

#ifdef ON_FULL_ILLUM_ARRAYS
  mem_protected::FreeMemHost(grid_host.energy);
  mem_protected::FreeMemHost(grid_host.stream);
  mem_protected::FreeMemHost(grid_host.impuls);
#endif

  WRITE_LOG("Free host arrays\n");
}

#endif //! SEPARATE_GPU
