#ifdef USE_CUDA

#include "cuda_init_mem.h"
#include "cuda_memory.h"

cuda::geo::device_host_ptr_t device_host_ptr; //связь хоста с массивами внутри видеокарты (через эту структуры идёт обращение к массивам на карте)

namespace mem = cuda::mem_protected;

void cuda::InitDirectionsOnDevice(const grid_directions_t &grid_host, geo::grid_directions_device_t *&grid_device) {

  const int n = grid_host.size;

  DIE_IF(sizeof(geo::direction_device_t) != sizeof(grid_host.directions[0])); //размеры структур на хосте и видеокарте должны быть соразмерны

  mem::Malloc(sizeof(geo::grid_directions_device_t), &grid_device);              // память под структуру сетки направлений
  mem::Malloc(n * sizeof(geo::direction_device_t), &device_host_ptr.directions); //память на массив направлений внутри структуры сетки

  mem::CpyToDevice(&grid_device->size, &n, sizeof(n));                                                            // копируем размер сетки
  mem::CpyToDevice(&grid_device->full_area, &grid_host.full_area, sizeof(grid_host.full_area));                   //копируем площадь
  mem::CpyToDevice(&grid_device->directions, &device_host_ptr.directions, sizeof(device_host_ptr.directions));    //указатель со стороны хоста в указатель на стороне карты
  mem::CpyToDevice(device_host_ptr.directions, grid_host.directions.data(), n * sizeof(grid_host.directions[0])); // через связанный с хостом указатель копируем данные в grid_device
}

void cuda::ClearDirectionsOnDevice(geo::grid_directions_device_t *&grid_device) {
  mem::FreeMem(device_host_ptr.directions); // удаляем массив внутри структуры
  mem::FreeMem(grid_device);                // удаляем саму структуру
}

void cuda::InitGridOnDevice(const grid_t &grid_host,
                            const int size_dir, const int start_dir, const int end_dir,
                            geo::grid_device_t *&grid_device) {

  const int size_grid = grid_host.size;
  const int num_cell_loc = grid_host.loc_size;
  const int cell_loc_shift = grid_host.loc_shift;
  const int loc_size_dir = end_dir - start_dir;

  mem::Malloc(sizeof(geo::grid_device_t), (void **)&grid_device); // память под структуру сетки

  mem::CpyToDevice(&grid_device->size, &size_grid, sizeof(size_grid));                        // копируем размер сетки
  mem::CpyToDevice(&grid_device->local_scattering_size, &loc_size_dir, sizeof(loc_size_dir)); // копируем локальное количество направлений
  mem::CpyToDevice(&grid_device->local_scattering_disp, &start_dir, sizeof(start_dir));       // копируем начальное направление

  mem::CpyToDevice(&grid_device->loc_size, &num_cell_loc, sizeof(num_cell_loc));  // копируем локальный размер сетки
  mem::CpyToDevice(&grid_device->shift, &cell_loc_shift, sizeof(cell_loc_shift)); // копируем начало локальной сетки

  //излучение
  mem::Malloc((CELL_SIZE * size_dir * size_grid * sizeof(Type)), &device_host_ptr.illum);       //память на массив излучения внутри структуры сетки
  mem::CpyToDevice(&grid_device->illum, &device_host_ptr.illum, sizeof(device_host_ptr.illum)); //указатель со стороны хоста в указатель на стороне карты

  //интеграл рассеяния
  mem::Malloc(loc_size_dir * size_grid * sizeof(Type), &device_host_ptr.int_scattering);                                   //память на массив рассеяния внутри структуры сетки
  mem::CpyToDevice(&grid_device->int_scattering, &device_host_ptr.int_scattering, sizeof(device_host_ptr.int_scattering)); //указатель со стороны хоста в указатель на стороне карты

#ifndef ONLY_CUDA_SCATTERING
  //дивергенция потока
  mem::Malloc(num_cell_loc * sizeof(Type), &device_host_ptr.divstream);                                     //память на массив дивергенций внутри структуры сетки
  mem::CpyToDevice(&grid_device->divstream, &device_host_ptr.divstream, sizeof(device_host_ptr.divstream)); //указатель со стороны хоста в указатель на стороне карты

  // дивергенция импульса
  mem::Malloc(num_cell_loc * sizeof(Vector3), &device_host_ptr.divimpuls);                                  //память на массив дивергенций внутри структуры сетки
  mem::CpyToDevice(&grid_device->divimpuls, &device_host_ptr.divimpuls, sizeof(device_host_ptr.divimpuls)); //указатель со стороны хоста в указатель на стороне карты

  //объемы
  mem::Malloc(size_grid * sizeof(Type), &device_host_ptr.volume);                                  //память на массив объемов внутри структуры сетки
  mem::CpyToDevice(&grid_device->volume, &device_host_ptr.volume, sizeof(device_host_ptr.volume)); //указатель со стороны хоста в указатель на стороне карты

  //площади
  mem::Malloc(CELL_SIZE * size_grid * sizeof(Type), &device_host_ptr.areas);                    //память на массив площадей внутри структуры сетки
  mem::CpyToDevice(&grid_device->areas, &device_host_ptr.areas, sizeof(device_host_ptr.areas)); //указатель со стороны хоста в указатель на стороне карты

  //нормали
  mem::Malloc(CELL_SIZE * size_grid * sizeof(Vector3), &device_host_ptr.normals);                     //память на массив нормалей внутри структуры сетки
  mem::CpyToDevice(&grid_device->normals, &device_host_ptr.normals, sizeof(device_host_ptr.normals)); //указатель со стороны хоста в указатель на стороне карты

#ifdef ON_FULL_ILLUM_ARRAYS
  // энергия
  mem::Malloc(num_cell_loc * sizeof(Type), &device_host_ptr.energy);                               //память на массив энергии внутри структуры сетки
  mem::CpyToDevice(&grid_device->energy, &device_host_ptr.energy, sizeof(device_host_ptr.energy)); //указатель со стороны хоста в указатель на стороне карты

  //поток
  mem::Malloc(num_cell_loc * sizeof(Vector3), &device_host_ptr.stream);                            //память на массив потока внутри структуры сетки
  mem::CpyToDevice(&grid_device->stream, &device_host_ptr.stream, sizeof(device_host_ptr.stream)); //указатель со стороны хоста в указатель на стороне карты

  //импульс
  mem::Malloc(num_cell_loc * sizeof(Matrix3), &device_host_ptr.impuls);                            //память на массив импульса внутри структуры сетки
  mem::CpyToDevice(&grid_device->impuls, &device_host_ptr.impuls, sizeof(device_host_ptr.impuls)); //указатель со стороны хоста в указатель на стороне карты

#endif
#endif
}

void cuda::ClearGridOnDevice(geo::grid_device_t *&grid_device) {

  mem::FreeMem(device_host_ptr.illum);
  mem::FreeMem(device_host_ptr.int_scattering);
  mem::FreeMem(device_host_ptr.divimpuls);
  mem::FreeMem(device_host_ptr.divstream);

  mem::FreeMem(device_host_ptr.normals);
  mem::FreeMem(device_host_ptr.volume);
  mem::FreeMem(device_host_ptr.areas);

#ifdef ON_FULL_ILLUM_ARRAYS
  mem::FreeMem(device_host_ptr.energy);
  mem::FreeMem(device_host_ptr.stream);
  mem::FreeMem(device_host_ptr.impuls);
#endif

  mem::FreeMem(grid_device);
}

#endif //! USE_CUDA