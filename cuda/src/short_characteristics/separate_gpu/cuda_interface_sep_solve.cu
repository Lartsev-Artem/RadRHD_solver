#ifdef USE_CUDA

#include "cuda_def.h"
#include "cuda_illum_sep_param.h"
#include "cuda_init_mem.h"
#include "cuda_memory.h"

#include "cuda_interface.h"
#include "cuda_multi_interface.h"

#include "cuda_illum_param.h"
#include "cuda_scattering.h"

#include "cuda_struct.h"

using namespace cuda::geo;

#ifdef SEPARATE_GPU

#ifdef MULTI_GPU
int cuda::interface::separate_device::CalculateIntScatteringAsync(const grid_directions_t &grid_dir, grid_t &grid,
                                                                  const IdType start_dir, const IdType end_dir, const e_cuda_stream_id_t stream) {
  const IdType M_loc = end_dir - start_dir;
  const IdType M = grid_dir.size;

  for (int id_dev = 0; id_dev < gpu_config.GPU_N; id_dev++) {

    const IdType N = gpu_config.size[id_dev];
    const IdType dispN = gpu_config.disp[id_dev];
    CUDA_TREADS_2D(threads);
    CUDA_BLOCKS_2D(blocks, N, M_loc);

    CUDA_CALL_FUNC(cudaSetDevice, id_dev);

    cudaMemcpyAsync(device_host_ptrN[id_dev].illum, grid.Illum + M * dispN, M * N * sizeof(grid.Illum[0]), cudaMemcpyHostToDevice, cuda_streams[stream]);

    kernel::GetS_MPI_multi_device<<<blocks, threads, 0, cuda_streams[stream]>>>(grid_dir_deviceN[id_dev], grid_deviceN[id_dev], N, start_dir, end_dir);

    CUDA_CALL_FUNC(cudaMemcpyAsync, grid.scattering + dispN, device_host_ptrN[id_dev].int_scattering, M_loc * N * sizeof(grid.scattering[0]), cudaMemcpyDeviceToHost, cuda_streams[stream]);
  }
}
#elif defined SINGLE_GPU
int cuda::interface::separate_device::CalculateIntScatteringAsync(const grid_directions_t &grid_dir, grid_t &grid,
                                                                  const IdType start_dir, const IdType end_dir, const e_cuda_stream_id_t stream) {

  const IdType M_loc = end_dir - start_dir;
  const IdType M = grid_dir.size;
  constexpr int id_dev = 0;

  for (int it = 0; it < GPU_DIV_PARAM; it++) {

    const IdType N = gpu_config.size[it];
    const IdType dispN = gpu_config.disp[it];
    CUDA_TREADS_2D(threads);
    CUDA_BLOCKS_2D(blocks, N, M_loc);

    cudaMemcpyAsync(device_host_ptrN[id_dev].illum, grid.Illum + M * dispN, M * N * sizeof(grid.Illum[0]),
                    cudaMemcpyHostToDevice, cuda_streams[stream]);

    CudaSyncStream(e_cuda_scattering_2);

    kernel::GetS_MPI_multi_device<<<blocks, threads, 0, cuda_streams[stream]>>>(grid_dir_deviceN[id_dev], grid_deviceN[id_dev], N, start_dir, end_dir);

    CUDA_CALL_FUNC(cudaGetLastError);

    IdType size = M_loc * N * sizeof(grid.scattering[0]);

    CudaSyncStream(stream);

    CUDA_CALL_FUNC(cudaMemcpyAsync, grid.scattering + M_loc * dispN, device_host_ptrN[id_dev].int_scattering, size, cudaMemcpyDeviceToHost, cuda_streams[e_cuda_scattering_2]);

    {
      CudaSyncStream(e_cuda_params);
      CudaSyncStream(e_cuda_scattering_2);

      /// \todo: правильный признак
      if (grid.energy != nullptr) {
        cuda::interface::separate_device::CalculateAllParamAsync(0, it, grid_dir, grid, e_cuda_params);
      }
      CudaSyncStream(e_cuda_params);
      CudaSyncStream(e_cuda_scattering_2);
    }
  }

  return e_completion_success;
}
#endif

#if 0 // DBG
int cuda::interface::separate_device::CalculateIntScatteringAsync(const grid_directions_t &grid_dir, grid_t &grid,
                                                                  const IdType start_dir, const IdType end_dir, const e_cuda_stream_id_t stream) {
  const IdType M = end_dir - start_dir;
  const IdType N = grid.size;

  const int id_dev = 0;

  CUDA_TREADS_2D(threads);
  CUDA_BLOCKS_2D(blocks, N, M);

  cudaMemcpyAsync(device_host_ptrN[id_dev].illum, grid.Illum, grid_dir.size * N * sizeof(grid.Illum[0]),
                  cudaMemcpyHostToDevice, cuda_streams[stream]);

  kernel::GetS_MPI_multi_device<<<blocks, threads, 0, cuda_streams[stream]>>>(grid_dir_deviceN[id_dev], grid_deviceN[id_dev], grid.size - 0, start_dir, end_dir);

  CUDA_CALL_FUNC(cudaGetLastError);

  IdType disp = start_dir * N;
  IdType size = (end_dir - start_dir) * N * sizeof(grid.scattering[0]);

  CUDA_CALL_FUNC(cudaMemcpyAsync, grid.scattering + disp, device_host_ptrN[id_dev].int_scattering + disp, size, cudaMemcpyDeviceToHost, cuda_streams[stream]);
  return e_completion_success;
}
#endif

int cuda::interface::separate_device::CalculateAllParamAsync(const int id_dev, const int im_dev, const grid_directions_t &grid_dir, grid_t &grid, e_cuda_stream_id_t st) {

#ifdef ON_FULL_ILLUM_ARRAYS
#ifndef ONLY_CUDA_SCATTERING
  const IdType N_loc = gpu_config.size_params[im_dev];
  const IdType host_disp = gpu_config.disp_params[im_dev];
  // WRITE_LOG("N_loc=%ld,host_disp=%ld\n", N_loc, host_disp);

  CUDA_TREADS_1D(threads);
  CUDA_BLOCKS_1D(blocks, N_loc);

  cuda::separate_device::kernel::MakeIllumParam<<<blocks, threads, 0, cuda_streams[st]>>>(grid_dir_deviceN[id_dev], grid_deviceN[id_dev], N_loc, host_disp);

  CUDA_CALL_FUNC(cudaGetLastError);

  CUDA_CALL_FUNC(cudaMemcpyAsync, grid.energy + host_disp, device_host_ptrN[id_dev].energy, N_loc * sizeof(grid.energy[0]), cudaMemcpyDeviceToHost, cuda_streams[st]);
  CUDA_CALL_FUNC(cudaMemcpyAsync, grid.stream + host_disp, device_host_ptrN[id_dev].stream, N_loc * sizeof(grid.stream[0]), cudaMemcpyDeviceToHost, cuda_streams[st]);
  CUDA_CALL_FUNC(cudaMemcpyAsync, grid.impuls + host_disp, device_host_ptrN[id_dev].impuls, N_loc * sizeof(grid.impuls[0]), cudaMemcpyDeviceToHost, cuda_streams[st]);
#endif
#endif

  return e_completion_success;
}

#ifdef MULTI_GPU
int cuda::interface::separate_device::CalculateAllParamAsync(const grid_directions_t &grid_dir, grid_t &grid, e_cuda_stream_id_t st) {

  for (int id_dev = 0; id_dev < GPU_DIV_PARAM; id_dev++) {

    CUDA_CALL_FUNC(cudaSetDevice, id_dev);

    CudaSyncStream(e_cuda_scattering_1); //ждем излучения
    CudaSyncStream(st);                  //теоретически ждем предыдущую итерацию

    CalculateAllParamAsync(id_dev, id_dev, grid_dir, grid, st);
  }

  return e_completion_success;
}
#endif

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