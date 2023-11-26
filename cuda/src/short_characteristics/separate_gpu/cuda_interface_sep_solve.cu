#ifdef USE_CUDA

#include "cuda_def.h"
#include "cuda_init_mem.h"
#include "cuda_memory.h"

#include "cuda_interface.h"

#include "cuda_illum_param.h"
#include "cuda_scattering.h"

#include "cuda_struct.h"

using namespace cuda::geo;

#ifdef SEPARATE_GPU
int cuda::interface::CalculateIntScatteringAsync(const grid_directions_t &grid_dir, grid_t &grid,
                                                 const IdType start_dir, const IdType end_dir, const e_cuda_stream_id_t stream) {
  const IdType M = end_dir - start_dir;
  const IdType N = grid.size;

  CUDA_TREADS_2D(threads);
  CUDA_BLOCKS_2D(blocks, N, M);

  cudaMemcpyAsync(device_host_ptr.illum, grid.Illum, N * sizeof(grid.Illum[0]),
                  cudaMemcpyHostToDevice, cuda_streams[e_cuda_scattering_1]);

  kernel::GetS_MPI_multi_device<<<blocks, threads, 0, cuda_streams[stream]>>>(grid_dir_device, grid_device, grid.size - 0, start_dir, end_dir);

  CUDA_CALL_FUNC(cudaGetLastError);

  IdType disp = start_dir * N;
  IdType size = (end_dir - start_dir) * N * sizeof(grid.scattering[0]);

  CUDA_CALL_FUNC(cudaMemcpyAsync, grid.scattering + disp, device_host_ptr.int_scattering + disp, size, cudaMemcpyDeviceToHost, cuda_streams[stream]);

  return e_completion_success;
}

int cuda::interface::CalculateAllParamAsync(const grid_directions_t &grid_dir, grid_t &grid, e_cuda_stream_id_t st) {
  /// \todo реализация
  return 0;
}

#endif //! SEPARATE_GPU

#endif //! USE_CUDA

#if 0
void cuda::interface::CudaSendIllumAsyncMultiDev(const IdType size, const IdType shift, const Type *Illum_host) {

  D_LD;

  int Ndev = 3; ///< число карт
  CUDA_CALL_FUNC(cudaGetDeviceCount, &Ndev);

  int disp[Ndev] = {0, 1, 2};
  int send[Ndev] = {0, 1, 2};

  int M = 100; // dirs

  for (size_t i = 0; i < Ndev; i++) {

    cudaSetDevice(i);

    cudaMemcpyAsync(device_host_ptr.illum, Illum_host + disp[i], send[i] * sizeof(Illum_host[0]),
                    cudaMemcpyHostToDevice, cuda_streams[e_cuda_scattering_1]);

    CUDA_TREADS_2D(threads);
    CUDA_BLOCKS_2D(blocks, send[i], M);

    // kernel::GetS_MPI_multy_device<<<blocks, threads>>>(grid_dir_device, grid_device);

    // CUDA_CALL_FUNC(cudaGetLastError);

    cuda::mem_protected::CpyToHostAsync(Illum_host + disp[i], device_host_ptr.illum, send[i] * sizeof(Illum_host[0]));

    // sync if dev=1 and split config

}

for (int i = 0; i < Ndev; i++) {

  // Set device
  CUDA_CALL_FUNC(cudaSetDevice, i);

  // Wait for all operations to finish
  // cudaStreamSynchronize(plan[i].stream);
}
}
#endif