#ifdef USE_CUDA
#include "cuda_interface.h"

#include "cuda_init_mem.h"

#include "global_value.h"
#include "reader_bin.h"

#include "cuda_memory.h"

cuda::geo::grid_directions_device_t *grid_dir_device;
cuda::geo::grid_device_t *grid_device;

int cuda::interface::InitDevice(const std::string &address, const grid_directions_t &grid_dir_host, grid_t &grid_host) {

  CUDA_CALL_FUNC(cudaSetDevice, 0);

  InitDirectionsOnDevice(grid_dir_host, grid_dir_device);
  InitGridOnDevice(grid_host, grid_dir_host.size, grid_dir_host.loc_shift, grid_dir_host.loc_shift + grid_dir_host.loc_size, grid_device); /// \todo сдвиги в сетку и сюда сразу передавать геометрию

  /// \todo это отдельно, т.к. относится к инициализации хоста
  mem_protected::MallocHost((CELL_SIZE * grid_dir_host.size * grid_host.size * sizeof(Type)), &grid_host.Illum);
  mem_protected::MallocHost(((grid_dir_host.loc_size) * grid_host.size * sizeof(Type)), &grid_host.scattering);

  if (grid_dir_host.loc_shift == 0 || (grid_host.size != grid_host.loc_size)) // или нулевой узел, или данные разделены
  {
    mem_protected::MallocHost((grid_host.loc_size * sizeof(Type)), &grid_host.divstream);
    mem_protected::MallocHost((grid_host.loc_size * sizeof(Vector3)), &grid_host.divimpuls);

#ifdef ON_FULL_ILLUM_ARRAYS
    mem_protected::MallocHost((grid_host.loc_size * sizeof(Type)), &grid_host.energy);
    mem_protected::MallocHost((grid_host.loc_size * sizeof(Vector3)), &grid_host.stream);
    mem_protected::MallocHost((grid_host.loc_size * sizeof(Matrix3)), &grid_host.impuls);
#endif

    const int N = grid_host.size;
    std::vector<Normals> normals;
    std::vector<Type> areas_faces;
    std::vector<Type> volume;

    uint32_t err = 0;
    err |= files_sys::bin::ReadSimple(address + F_AREAS, areas_faces);
    err |= files_sys::bin::ReadNormals(address + F_NORMALS, normals);
    err |= files_sys::bin::ReadSimple(address + F_VOLUME, volume);
    if (err)
      return e_completion_fail;

    std::vector<Vector3> dev_norm(normals.size() * CELL_SIZE);
    for (size_t i = 0; i < normals.size(); i++) {
      for (size_t j = 0; j < CELL_SIZE; j++) {
        dev_norm[i * CELL_SIZE + j] = normals[i].n[j];
      }
    }

    mem_protected::CpyToDevice(device_host_ptr.normals, dev_norm.data(), dev_norm.size() * sizeof(dev_norm[0]));
    mem_protected::CpyToDevice(device_host_ptr.volume, volume.data(), N * sizeof(volume[0]));
    mem_protected::CpyToDevice(device_host_ptr.areas, areas_faces.data(), CELL_SIZE * N * sizeof(areas_faces[0]));
  }

  return e_completion_success;
}

void cuda::interface::ClearDevice() {

  ClearDirectionsOnDevice(grid_dir_device);
  ClearGridOnDevice(grid_device);

  // CUDA_CALL_FUNC(cudaDeviceReset);

  WRITE_LOG("Free device arrays\n");
}
void cuda::interface::ClearHost(grid_t &grid_host) {

  mem_protected::FreeMemHost(grid_host.Illum);
  mem_protected::FreeMemHost(grid_host.scattering);

  mem_protected::FreeMemHost(grid_host.divstream);
  mem_protected::FreeMemHost(grid_host.divimpuls);

#ifdef ON_FULL_ILLUM_ARRAYS
  mem_protected::FreeMemHost(grid_host.energy);
  mem_protected::FreeMemHost(grid_host.stream);
  mem_protected::FreeMemHost(grid_host.impuls);
#endif

  WRITE_LOG("Free host arrays\n");
}

#endif //! USE_CUDA