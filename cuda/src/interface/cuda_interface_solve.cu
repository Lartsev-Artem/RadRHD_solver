#ifdef USE_CUDA

#include "cuda_def.h"
#include "cuda_init_mem.h"
#include "cuda_memory.h"

#include "cuda_interface.h"

#include "cuda_illum_param.h"
#include "cuda_scattering.h"

#include "cuda_struct.h"

using namespace cuda::geo;

int cuda::interface::CalculateAllParam(const grid_directions_t &grid_dir, grid_t &grid) {

  const int N = grid.size;
  const int N_loc = grid.loc_size;

  CUDA_TREADS_1D(threads);
  CUDA_BLOCKS_1D(blocks, N);

  kernel::MakeIllumParam<<<blocks, threads>>>(grid_dir_device, grid_device);

  CUDA_CALL_FUNC(cudaGetLastError);

  // не асинхронно т.к. сразу затем идет расчет газовый
  mem_protected::CpyToHost(grid.divstream, device_host_ptr.divstream, N_loc * sizeof(grid.divstream[0]));
  mem_protected::CpyToHost(grid.divimpuls, device_host_ptr.divimpuls, N_loc * sizeof(grid.divimpuls[0]));

#ifdef ON_FULL_ILLUM_ARRAYS
  mem_protected::CpyToHostAsync(grid.energy, device_host_ptr.energy, N_loc * sizeof(grid.energy[0]));
  mem_protected::CpyToHostAsync(grid.stream, device_host_ptr.stream, N_loc * sizeof(grid.stream[0]));
  mem_protected::CpyToHostAsync(grid.impuls, device_host_ptr.impuls, N_loc * sizeof(grid.impuls[0]));
#endif

  return e_completion_success;
}

int cuda::interface::CalculateIntScattering(const grid_directions_t &grid_dir, grid_t &grid) {
  const int M = grid_dir.size;
  const int N = grid.size;

  mem_protected::CpyToDevice(device_host_ptr.illum, grid.Illum, N * M * CELL_SIZE * sizeof(grid.Illum[0]));

  CUDA_TREADS_2D(threads);
  CUDA_BLOCKS_2D(blocks, N, M);

  kernel::GetS<<<blocks, threads>>>(grid_dir_device, grid_device);

  CUDA_CALL_FUNC(cudaGetLastError);

  mem_protected::CpyToHost(grid.scattering, device_host_ptr.int_scattering, N * M * sizeof(grid.scattering[0]));

  return e_completion_success;
}

#endif //! USE_CUDA