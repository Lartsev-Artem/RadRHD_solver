#if defined USE_CUDA && defined SPECTRUM

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

void cuda::interface::separate_device::SendVelocity(const grid_t &grid) {

  std::vector<Vector3> vel(grid.size);
  for (size_t i = 0; i < grid.size; i++) {
    vel[i] = grid.cells[i].phys_val.v;
  }

#ifdef SINGLE_GPU
  mem_protected::CpyToDevice(device_host_ptrN[0].velocity, vel.data(), grid.size * sizeof(vel[0]));
#else

  for (int it = 0; it < device_host_ptrN.size(); it++) {
    int disp = gpu_config.disp[it];
    int size = gpu_config.size[it];
    mem_protected::CpyToDevice(device_host_ptrN[it].velocity, vel.data() + disp, size * sizeof(vel[0]));
  }
#endif
}

#ifdef MULTI_GPU
int cuda::interface::separate_device::CalculateSpectrumIntScattering(const grid_directions_t &grid_dir, grid_t &grid,
                                                                     const IdType start_dir, const IdType end_dir, const e_cuda_stream_id_t stream) {

#error "There is no implementation"
}
#elif defined SINGLE_GPU
int cuda::interface::separate_device::CalculateSpectrumIntScattering(const grid_directions_t &grid_dir, grid_t &grid,
                                                                     const IdType start_dir, const IdType end_dir, const e_cuda_stream_id_t stream) {

  const IdType M_loc = end_dir - start_dir;
  const IdType M = grid_dir.size;
  constexpr int id_dev = 0;
  const Type frq = 0.5 * (grid.frq_grid[grid.cur_frq_id] + grid.frq_grid[grid.cur_frq_id + 1]);

  for (int it = 0; it < GPU_DIV_PARAM; it++) {

    cuda::separate_device::kernel::SetImDevice<<<1, 1>>>(grid_deviceN[id_dev],
                                                         gpu_config.size[it],
                                                         gpu_config.disp[it],
                                                         gpu_config.size_params[it],
                                                         gpu_config.disp_params[it]);

    CUDA_CALL_FUNC(cudaDeviceSynchronize);

    const IdType N = gpu_config.size[it];
    const IdType dispN = gpu_config.disp[it];
    CUDA_TREADS_2D(threads);
    CUDA_BLOCKS_2D(blocks, N, M_loc);

    cudaMemcpyAsync(device_host_ptrN[id_dev].illum, grid.Illum + M * dispN, M * N * sizeof(grid.Illum[0]),
                    cudaMemcpyHostToDevice, cuda_streams[stream]);

    CudaSyncStream(e_cuda_scattering_2);

    kernel::Get_spectrum_multi_device<<<blocks, threads, 0, cuda_streams[stream]>>>(grid_dir_deviceN[id_dev], grid_deviceN[id_dev],
                                                                                    frq, end_dir);

    // kernel::GetS_MPI_multi_device<<<blocks, threads, 0, cuda_streams[stream]>>>(grid_dir_deviceN[id_dev], grid_deviceN[id_dev], N, start_dir, end_dir);

    CUDA_CALL_FUNC(cudaGetLastError);

    IdType size = M_loc * N * sizeof(grid.scattering[0]);

    CudaSyncStream(stream);

    CUDA_CALL_FUNC(cudaMemcpyAsync, grid.scattering + M_loc * dispN, device_host_ptrN[id_dev].int_scattering, size, cudaMemcpyDeviceToHost, cuda_streams[e_cuda_scattering_2]);

    CudaSyncStream(e_cuda_params);
    CudaSyncStream(e_cuda_scattering_2);
  }

  return e_completion_success;
}

int cuda::interface::separate_device::CalculateFullSpectrumIntScattering(const grid_directions_t &grid_dir, grid_t &grid,
                                                                         const IdType start_dir, const IdType end_dir, const e_cuda_stream_id_t stream) {

  const IdType M_loc = end_dir - start_dir;
  const IdType M = grid_dir.size;
  const IdType F = grid.size_frq;
  constexpr int id_dev = 0;

  for (int it = 0; it < GPU_DIV_PARAM; it++) {

    cuda::separate_device::kernel::SetImDevice<<<1, 1>>>(grid_deviceN[id_dev],
                                                         gpu_config.size[it],
                                                         gpu_config.disp[it],
                                                         gpu_config.size_params[it],
                                                         gpu_config.disp_params[it]);

    CUDA_CALL_FUNC(cudaDeviceSynchronize);

    const IdType N = gpu_config.size[it];
    const IdType dispN = gpu_config.disp[it];
    CUDA_TREADS_3D(threads);
    CUDA_BLOCKS_3D(blocks, N, M_loc, F);

    cudaMemcpyAsync(device_host_ptrN[id_dev].illum, grid.Illum + M * dispN * F, M * N * F * sizeof(grid.Illum[0]),
                    cudaMemcpyHostToDevice, cuda_streams[stream]);

    CudaSyncStream(e_cuda_scattering_2);

    kernel::Get_full_spectrum_multi_device<<<blocks, threads, 0, cuda_streams[stream]>>>(grid_dir_deviceN[id_dev], grid_deviceN[id_dev],
                                                                                         end_dir);

    // kernel::GetS_MPI_multi_device<<<blocks, threads, 0, cuda_streams[stream]>>>(grid_dir_deviceN[id_dev], grid_deviceN[id_dev], N, start_dir, end_dir);

    CUDA_CALL_FUNC(cudaGetLastError);

    IdType size = F * M_loc * N * sizeof(grid.scattering[0]);

    CudaSyncStream(stream);

    CUDA_CALL_FUNC(cudaMemcpyAsync, grid.scattering + M_loc * dispN * F, device_host_ptrN[id_dev].int_scattering, size, cudaMemcpyDeviceToHost, cuda_streams[e_cuda_scattering_2]);

    CudaSyncStream(e_cuda_params);
    CudaSyncStream(e_cuda_scattering_2);
  }

  return e_completion_success;
}
#endif
#endif
#endif