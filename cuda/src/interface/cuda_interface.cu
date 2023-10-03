#ifdef USE_CUDA
#include "cuda_interface.h"

#include "cuda_init_mem.h"
#include "cuda_struct.h"

#include "global_value.h"
#include "reader_bin.h"

#include "cuda_memory.h"

cuda::geo::grid_directions_device_t *grid_dir_device;
cuda::geo::grid_device_t *grid_device;

int InitDevice(const std::string &address, const grid_directions_t &grid_dir_host, grid_t &grid_host, const int start, const int end) {

  cuda::InitDirectionsOnDevice(grid_dir_host, grid_dir_device);
  cuda::InitGridOnDevice(grid_host, grid_dir_host.size, start, end, grid_device); /// \todo сдвиги в сетку и сюда сразу передавать геометрию

  /// \todo это отдельно, т.к. относится к инициализации хоста
  cuda::mem_protected::MallocHost((CELL_SIZE * grid_dir_host.size * grid_host.size * sizeof(Type)), &grid_host.Illum);
  cuda::mem_protected::MallocHost(((end - start) * grid_host.size * sizeof(Type)), &grid_host.scattering);

  if (start == 0 || (grid_host.size != grid_host.loc_size)) // или нулевой узел, или данные разделены
  {
    cuda::mem_protected::MallocHost((grid_host.loc_size * sizeof(Type)), &grid_host.divstream);
    cuda::mem_protected::MallocHost((grid_host.loc_size * sizeof(Vector3)), &grid_host.divimpuls);

#ifdef ON_FULL_ILLUM_ARRAYS
    cuda::mem_protected::MallocHost((grid_host.loc_size * sizeof(Type)), &grid_host.energy);
    cuda::mem_protected::MallocHost((grid_host.loc_size * sizeof(Vector3)), &grid_host.stream);
    cuda::mem_protected::MallocHost((grid_host.loc_size * sizeof(Matrix3)), &grid_host.impuls);
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

    cuda::mem_protected::CpyToDevice(cuda::device_host_ptr.normals, dev_norm.data(), dev_norm.size() * sizeof(dev_norm[0]));
    cuda::mem_protected::CpyToDevice(cuda::device_host_ptr.volume, volume.data(), N * sizeof(volume[0]));
    cuda::mem_protected::CpyToDevice(cuda::device_host_ptr.areas, areas_faces.data(), CELL_SIZE * N * sizeof(areas_faces[0]));
  }

  return e_completion_success;
}
#endif //! USE_CUDA