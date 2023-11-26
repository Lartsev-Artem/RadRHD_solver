#ifdef USE_CUDA

#include "cuda_def.h"
#include "cuda_init_mem.h"
#include "cuda_memory.h"

#include "cuda_interface.h"

#include "cuda_illum_param.h"
#include "cuda_scattering.h"

#include "cuda_struct.h"

using namespace cuda::geo;

#ifndef SEPARATE_GPU

int cuda::interface::CalculateAllParam(const grid_directions_t &grid_dir, grid_t &grid) {

#ifndef ONLY_CUDA_SCATTERING

  const int N_loc = grid.loc_size;

  CUDA_TREADS_1D(threads);
  CUDA_BLOCKS_1D(blocks, N_loc);

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
#endif
  return e_completion_success;
}

int cuda::interface::CalculateIntScattering(const grid_directions_t &grid_dir, grid_t &grid) {
  const IdType M = grid_dir.size;
  const IdType N = grid.size;

  mem_protected::CpyToDevice(device_host_ptr.illum, grid.Illum, N * M * CELL_SIZE * sizeof(grid.Illum[0]));

  CUDA_TREADS_2D(threads);
  CUDA_BLOCKS_2D(blocks, N, M);

  kernel::GetS<<<blocks, threads>>>(grid_dir_device, grid_device);

  CUDA_CALL_FUNC(cudaGetLastError);

  mem_protected::CpyToHost(grid.scattering, device_host_ptr.int_scattering, N * M * sizeof(grid.scattering[0]));

  return e_completion_success;
}

int cuda::interface::CalculateIntScatteringAsync(const grid_directions_t &grid_dir, grid_t &grid,
                                                 const IdType start_dir, const IdType end_dir, const e_cuda_stream_id_t stream) {
  const IdType M = end_dir - start_dir; // grid_dir.size;
  const IdType N = grid.size;

  CUDA_TREADS_2D(threads);
  CUDA_BLOCKS_2D(blocks, N, M);

  // dim3 threads(32, 16);
  // dim3 blocks((N + 32 - 1) / 32, (M + 16 - 1) / 16);

  // надо как то дать задержку второму потоку, относительно первого
  kernel::GetS_MPI_Stream<<<blocks, threads, 0, cuda_streams[stream]>>>(grid_dir_device, grid_device, start_dir, end_dir);

  CUDA_CALL_FUNC(cudaGetLastError);

  IdType disp = start_dir * N;
  IdType size = (end_dir - start_dir) * N * sizeof(grid.scattering[0]);

  CUDA_CALL_FUNC(cudaMemcpyAsync, grid.scattering + disp, device_host_ptr.int_scattering + disp, size, cudaMemcpyDeviceToHost, cuda_streams[stream]);

  return e_completion_success;
}
int cuda::interface::CalculateAllParamAsync(const grid_directions_t &grid_dir, grid_t &grid, e_cuda_stream_id_t st) {

#ifndef ONLY_CUDA_SCATTERING
  const IdType N_loc = grid.loc_size;

  CUDA_TREADS_1D(threads);
  CUDA_BLOCKS_1D(blocks, N_loc);

  kernel::MakeIllumParam<<<blocks, threads, 0, cuda_streams[st]>>>(grid_dir_device, grid_device);

  CUDA_CALL_FUNC(cudaGetLastError);

  CUDA_CALL_FUNC(cudaMemcpyAsync, grid.divstream, device_host_ptr.divstream, N_loc * sizeof(grid.divstream[0]), cudaMemcpyDeviceToHost, cuda_streams[st]);
  CUDA_CALL_FUNC(cudaMemcpyAsync, grid.divimpuls, device_host_ptr.divimpuls, N_loc * sizeof(grid.divimpuls[0]), cudaMemcpyDeviceToHost, cuda_streams[st]);

#ifdef ON_FULL_ILLUM_ARRAYS
  CUDA_CALL_FUNC(cudaMemcpyAsync, grid.energy, device_host_ptr.energy, N_loc * sizeof(grid.energy[0]), cudaMemcpyDeviceToHost, cuda_streams[st]);
  CUDA_CALL_FUNC(cudaMemcpyAsync, grid.stream, device_host_ptr.stream, N_loc * sizeof(grid.stream[0]), cudaMemcpyDeviceToHost, cuda_streams[st]);
  CUDA_CALL_FUNC(cudaMemcpyAsync, grid.impuls, device_host_ptr.impuls, N_loc * sizeof(grid.impuls[0]), cudaMemcpyDeviceToHost, cuda_streams[st]);
#endif
#endif

  return e_completion_success;
}

#endif //! SEPARATE_GPU

#endif //! USE_CUDA
