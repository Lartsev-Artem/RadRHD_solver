#include "cuda_memory.h"
#include "cuda_multi_init.h"

#ifdef SEPARATE_GPU
namespace sep = cuda::separate_device;
namespace mem = cuda::mem_protected;

void sep::InitDirectionsOnMultiDevice(const grid_directions_t &grid_host,
                                      geo::device_host_ptr_t &device_host_ptr,
                                      geo::grid_directions_device_t *&grid_device) {

  const IdType n = grid_host.size;

  DIE_IF(sizeof(geo::direction_device_t) != sizeof(grid_host.directions[0])); //размеры структур на хосте и видеокарте должны быть соразмерны

  mem::Malloc(sizeof(geo::grid_directions_device_t), &grid_device);              // память под структуру сетки направлений
  mem::Malloc(n * sizeof(geo::direction_device_t), &device_host_ptr.directions); //память на массив направлений внутри структуры сетки

  mem::CpyToDevice(&grid_device->size, &n, sizeof(n));                                                            // копируем размер сетки
  mem::CpyToDevice(&grid_device->full_area, &grid_host.full_area, sizeof(grid_host.full_area));                   //копируем площадь
  mem::CpyToDevice(&grid_device->directions, &device_host_ptr.directions, sizeof(device_host_ptr.directions));    //указатель со стороны хоста в указатель на стороне карты
  mem::CpyToDevice(device_host_ptr.directions, grid_host.directions.data(), n * sizeof(grid_host.directions[0])); // через связанный с хостом указатель копируем данные в grid_device
}

void sep::InitMultiDeviceGrid(int id_dev, const multi_gpu_config_t &gpu_conf, const grid_t &grid_host,
                              const grid_directions_t &grid_dir_host,
                              geo::device_host_ptr_t &device_host_ptr,
                              geo::grid_device_t *&grid_device) {

  const IdType size_grid = grid_host.size;
  const IdType loc_shift_grid = grid_host.loc_shift;
  const IdType loc_size_grid = grid_host.loc_size;

  const IdType num_cell_loc = gpu_conf.size[id_dev];   // grid_host.loc_size;
  const IdType cell_loc_shift = gpu_conf.disp[id_dev]; // grid_host.loc_shift;

  const IdType loc_size_dir = grid_dir_host.loc_size;
  const IdType dir_shift = grid_dir_host.loc_shift;

  mem::Malloc(sizeof(geo::grid_device_t), (void **)&grid_device); // память под структуру сетки

  mem::CpyToDevice(&grid_device->size, &size_grid, sizeof(size_grid));                        // копируем размер сетки
  mem::CpyToDevice(&grid_device->local_scattering_size, &loc_size_dir, sizeof(loc_size_dir)); // копируем локальное количество направлений
  mem::CpyToDevice(&grid_device->local_scattering_disp, &dir_shift, sizeof(dir_shift));       // копируем начальное направление

  mem::CpyToDevice(&grid_device->loc_size_gpu, &num_cell_loc, sizeof(num_cell_loc));  // копируем локальный размер сетки
  mem::CpyToDevice(&grid_device->shift_gpu, &cell_loc_shift, sizeof(cell_loc_shift)); // копируем начало локальной сетки

  mem::CpyToDevice(&grid_device->loc_size_params, &loc_size_grid, sizeof(loc_size_grid)); // копируем локальный размер сетки
  mem::CpyToDevice(&grid_device->shift_params, &loc_shift_grid, sizeof(loc_shift_grid));  // копируем начало локальной сетки

  //излучение
  mem::Malloc((grid_dir_host.size * num_cell_loc * sizeof(Type)), &device_host_ptr.illum);      //память на массив излучения внутри структуры сетки
  mem::CpyToDevice(&grid_device->illum, &device_host_ptr.illum, sizeof(device_host_ptr.illum)); //указатель со стороны хоста в указатель на стороне карты

  //интеграл рассеяния
  mem::Malloc(loc_size_dir * num_cell_loc * sizeof(Type), &device_host_ptr.int_scattering);                                //память на массив рассеяния внутри структуры сетки
  mem::CpyToDevice(&grid_device->int_scattering, &device_host_ptr.int_scattering, sizeof(device_host_ptr.int_scattering)); //указатель со стороны хоста в указатель на стороне карты

#ifdef ON_FULL_ILLUM_ARRAYS

  int n = gpu_conf.size_params[id_dev];
#ifdef SINGLE_GPU
  n = *std::max_element(gpu_conf.size_params.begin(), gpu_conf.size_params.end());
#endif

  // энергия
  mem::Malloc(n * sizeof(Type), &device_host_ptr.energy);                                          //память на массив энергии внутри структуры сетки
  mem::CpyToDevice(&grid_device->energy, &device_host_ptr.energy, sizeof(device_host_ptr.energy)); //указатель со стороны хоста в указатель на стороне карты

  //поток
  mem::Malloc(n * sizeof(Vector3), &device_host_ptr.stream);                                       //память на массив потока внутри структуры сетки
  mem::CpyToDevice(&grid_device->stream, &device_host_ptr.stream, sizeof(device_host_ptr.stream)); //указатель со стороны хоста в указатель на стороне карты

  //импульс
  mem::Malloc(n * sizeof(Matrix3), &device_host_ptr.impuls);                                       //память на массив импульса внутри структуры сетки
  mem::CpyToDevice(&grid_device->impuls, &device_host_ptr.impuls, sizeof(device_host_ptr.impuls)); //указатель со стороны хоста в указатель на стороне карты

#endif
}

void sep::ClearDirectionsOnMultiDevice(const multi_gpu_config_t &gpu_conf, std::vector<geo::device_host_ptr_t> &device_host_ptrs,
                                       std::vector<geo::grid_directions_device_t *> &grid_devices) {

  DIE_IF(device_host_ptrs.size() != grid_devices.size())
  for (size_t i = 0; i < grid_devices.size(); i++) {
    mem::FreeMem(device_host_ptrs[i].directions); // удаляем массив внутри структуры
    mem::FreeMem(grid_devices[i]);                // удаляем саму структуру
  }
}

void sep::ClearGridOnMultiDevice(multi_gpu_config_t &gpu_conf, std::vector<geo::device_host_ptr_t> &device_host_ptrs,
                                 std::vector<geo::grid_device_t *> &grid_devices) {

  DIE_IF(device_host_ptrs.size() != grid_devices.size())

  for (size_t i = 0; i < grid_devices.size(); i++) {
    mem::FreeMem(device_host_ptrs[i].illum);
    mem::FreeMem(device_host_ptrs[i].int_scattering);

#if defined ON_FULL_ILLUM_ARRAYS
    mem::FreeMem(device_host_ptrs[i].energy);
    mem::FreeMem(device_host_ptrs[i].stream);
    mem::FreeMem(device_host_ptrs[i].impuls);
#endif

    mem::FreeMem(grid_devices[i]);
  }

  gpu_conf.size.clear();
  gpu_conf.disp.clear();
  gpu_conf.GPU_N = 0;
}

#endif //! SEPARATE_GPU