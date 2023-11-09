#include "cuda_integrator.h"
#include "cuda_scattering.h"
#include "global_def.h"

__global__ void cuda::kernel::GetS_MPI_multy_device(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid,
                                                    const IdType startCell, const IdType endCell, const IdType start, const IdType end) {
  const IdType N = grid->size;

  const IdType i = blockIdx.x * blockDim.x + threadIdx.x;
  const IdType k = blockIdx.y * blockDim.y + threadIdx.y;

  const IdType size_loc = endCell - startCell;

  if (i >= size_loc || k >= end || k < start)
    return;

  const IdType M = dir->size;

  const Vector3 &cur_dir = dir->directions[grid->local_scattering_disp + k].dir;
  const Type *Illum = grid->illum;
  const geo::direction_device_t *all_dir = dir->directions;

  Type scatter = 0;
  for (IdType num_direction = 0; num_direction < M; num_direction++) {
    Type I = Illum[i * M + num_direction];
    scatter += device::direction_integrator::Gamma(all_dir[num_direction].dir, cur_dir) * I * all_dir[num_direction].area;
  }

  grid->int_scattering[k * N + i] = scatter / dir->full_area;
}

#include "cuda_def.h"
#include "cuda_interface.h"

#include "cuda_memory.h"
#include "cuda_struct.h"

void cuda::interface::CudaSendIllumAsync(const IdType size, const IdType shift, const Type *Illum_host) {

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

struct multi_gpu_config_t {
  int GPU_N;
  std::vector<IdType> size;
  std::vector<IdType> disp;
};
multi_gpu_config_t gpu_conf;

namespace mem = cuda::mem_protected;
#include "cuda_init_mem.h"

/// \todo rename + namespace multicuda
using namespace cuda;
void InitGridOnDeviceMultiCuda(int id_dev, const grid_t &grid_host,
                               const IdType size_dir, const IdType start_dir, const IdType end_dir,
                               geo::grid_device_t *&grid_device) {

  const IdType size_grid = grid_host.size;
  const IdType num_cell_loc = gpu_conf.size[id_dev];   // grid_host.loc_size;
  const IdType cell_loc_shift = gpu_conf.disp[id_dev]; // grid_host.loc_shift;
  const IdType loc_size_dir = end_dir - start_dir;

  mem::Malloc(sizeof(geo::grid_device_t), (void **)&grid_device); // память под структуру сетки

  mem::CpyToDevice(&grid_device->size, &size_grid, sizeof(size_grid));                        // копируем размер сетки
  mem::CpyToDevice(&grid_device->local_scattering_size, &loc_size_dir, sizeof(loc_size_dir)); // копируем локальное количество направлений
  mem::CpyToDevice(&grid_device->local_scattering_disp, &start_dir, sizeof(start_dir));       // копируем начальное направление

  mem::CpyToDevice(&grid_device->loc_size, &num_cell_loc, sizeof(num_cell_loc));  // копируем локальный размер сетки
  mem::CpyToDevice(&grid_device->shift, &cell_loc_shift, sizeof(cell_loc_shift)); // копируем начало локальной сетки

  //излучение
  mem::Malloc((size_dir * num_cell_loc * sizeof(Type)), &device_host_ptr.illum);                //память на массив излучения внутри структуры сетки
  mem::CpyToDevice(&grid_device->illum, &device_host_ptr.illum, sizeof(device_host_ptr.illum)); //указатель со стороны хоста в указатель на стороне карты

  //интеграл рассеяния
  mem::Malloc(loc_size_dir * num_cell_loc * sizeof(Type), &device_host_ptr.int_scattering);                                //память на массив рассеяния внутри структуры сетки
  mem::CpyToDevice(&grid_device->int_scattering, &device_host_ptr.int_scattering, sizeof(device_host_ptr.int_scattering)); //указатель со стороны хоста в указатель на стороне карты
}

#include "mpi_shifts.h"

std::vector<cuda::geo::grid_directions_device_t *> grid_dir_device;
std::vector<cuda::geo::grid_device_t *> grid_device;

int cuda::interface::InitDeviceExtraSize(const std::string &address, const grid_directions_t &grid_dir_host, grid_t &grid_host) {

  CUDA_CALL_FUNC(cudaGetDeviceCount, &gpu_conf.GPU_N);
  GetSend(gpu_conf.GPU_N, grid_host.size, gpu_conf.size);
  GetDisp(gpu_conf.GPU_N, grid_host.size, gpu_conf.disp);

  grid_dir_device.resize(gpu_conf.GPU_N);
  grid_device.resize(gpu_conf.GPU_N);

  for (int dev_id = 0; dev_id < gpu_conf.GPU_N; dev_id++) {
    CUDA_CALL_FUNC(cudaSetDevice, dev_id);
    InitDirectionsOnDevice(grid_dir_host, grid_dir_device[dev_id]);
    InitGridOnDeviceMultiCuda(dev_id, grid_host, grid_dir_host.size, grid_dir_host.loc_shift, grid_dir_host.loc_shift + grid_dir_host.loc_size, grid_device[dev_id]); /// \todo сдвиги в сетку и сюда сразу передавать геометрию
  }

  /// \todo это отдельно, т.к. относится к инициализации хоста
  mem_protected::MallocHost((grid_dir_host.size * grid_host.size * sizeof(Type)), &grid_host.Illum);
  mem_protected::MallocHost(((grid_dir_host.loc_size) * grid_host.size * sizeof(Type)), &grid_host.scattering);

  if (grid_dir_host.loc_shift == 0 || (grid_host.size != grid_host.loc_size)) // или нулевой узел, или данные разделены
  {
    mem_protected::MallocHost((grid_host.loc_size * sizeof(Type)), &grid_host.energy);
    mem_protected::MallocHost((grid_host.loc_size * sizeof(Vector3)), &grid_host.stream);
    mem_protected::MallocHost((grid_host.loc_size * sizeof(Matrix3)), &grid_host.impuls);

    // //сопоставимо с хранением energy и stream и impuls
    // const Type abs_op = grid.cells[cell].illum_val.absorp_coef;
    // const Type scat_op = grid.cells[cell].illum_val.scat_coef;
    // const Type rho = grid.cells[cell].phys_val.d;
    // const Type prs = grid.cells[cell].phys_val.p;
    // const Vector3 &v = grid.cells[cell].phys_val.v;
  }
  return e_completion_success;
}

__device__ Type IntegrateByCell(const IdType num_cell, const geo::grid_directions_device_t *dir, const geo::grid_device_t *grid) {
  const IdType M = dir->size;
  const IdType N = grid->size;

  Type res = 0;
  for (IdType i = 0; i < M; i++) {
    IdType pos = CELL_SIZE * (N * i + num_cell);

    Type I = 0;
    for (IdType k = 0; k < CELL_SIZE; k++) {
      I += grid->illum[pos + k];
    }
    I /= CELL_SIZE;

    res += I * dir->directions[i].area;
  }

  return res / dir->full_area;
}
__device__ void MakeEnergy(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid) {
  const IdType N = grid->loc_size;
  const IdType shift = grid->shift;
  const IdType i = blockIdx.x * blockDim.x + threadIdx.x;
  const IdType M = dir->size;

  if (i >= N)
    return;

  Type sum = 0;
  for (IdType k = 0; k < M; k++) {
    sum += grid->illum[M * i + k] * dir->directions[k].area;
  }

  grid->energy[i] = sum / dir->full_area; // direction_integrator::IntegrateByCell(shift + i, dir, grid);
}