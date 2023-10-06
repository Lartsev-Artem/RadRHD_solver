#if !defined CUDA_STRUCT_H && defined USE_CUDA
#define CUDA_STRUCT_H

// workaround issue between gcc >= 4.7 and cuda 5.5 (совместимость с компиляторами (см. док. Eigen3: https://eigen.tuxfamily.org/dox/TopicCUDA.html))
#if (defined __GNUC__) && (__GNUC__ > 4 || __GNUC_MINOR__ >= 7)
#undef _GLIBCXX_ATOMIC_BUILTINS
#undef _GLIBCXX_USE_INT128
#endif
#include <Eigen/Dense> /// требует компиляции с флагом --expt-relaxed-constexpr

typedef Eigen::Vector3d Vector3;
typedef Eigen::Matrix3d Matrix3;
typedef double Type;

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "dbgdef.h"
#include "solvers_config.h"

#include "geo_types.h" /// \todo: Eigen подключается там!!!

namespace cuda {
namespace geo {

/**
 * @brief  структура направления на сфере направлений
 *
 */
struct direction_device_t {
  Vector3 dir; ///< направление излучения
  Type area;   ///< площадь ячейке на сфере направлений
};

/**
 * @brief Сетка сферы направлений (аналог direction_grid_cpu)
 *
 */
struct grid_directions_device_t {

  int size;                       ///< число направлений
  direction_device_t *directions; ///< массив направлений
  Type full_area;                 ///< площадь сферы направлений

  __host__ __device__ grid_directions_device_t() : size(0), directions(nullptr), full_area(0) {}
  __host__ __device__ grid_directions_device_t(const grid_directions_t &grid_host);
  __host__ __device__ ~grid_directions_device_t() {}
};

/**
 * @brief структура сетки (аналог grid_t)
 *
 */
struct grid_device_t {
  int size; ///< размер всей сетки

  int loc_size; ///< локальное число ячеек (mpi)
  int shift;    ///< локальный сдвиг по ячейкам(mpi)

  int local_scattering_size; ///< локальный размер по направлениям(потоки)
  int local_scattering_disp; ///< локальный сдвиг  по направлениям(потоки)

  Type *illum;          ///< излучение
  Type *int_scattering; ///< интеграл рессеяния

  Type *divstream;    ///< дивергенция потока
  Vector3 *divimpuls; ///< дивергенция импульса

  Vector3 *normals; ///< нормали к ячейкам
  Type *areas;      ///< площадь граней
  Type *volume;     ///< объемы ячеек

#ifdef ON_FULL_ILLUM_ARRAYS
  Type *energy;    ///< энергия излучения
  Vector3 *stream; ///< поток энергии излучения
  Matrix3 *impuls; ///< импулься энергии излучения
#endif

  __host__ __device__ grid_device_t() : size(0), illum(nullptr), int_scattering(nullptr), divstream(nullptr),
                                        divimpuls(nullptr), normals(nullptr), areas(nullptr), volume(nullptr)
#ifdef ON_FULL_ILLUM_ARRAYS
                                        ,
                                        energy(nullptr), stream(nullptr), impuls(nullptr) {
  }
#endif
};

struct device_host_ptr_t {
  direction_device_t *directions;

  Type *illum;
  Type *int_scattering;

  Type *divstream;
  Vector3 *divimpuls;

  Vector3 *normals;
  Type *areas;
  Type *volume;

#ifdef ON_FULL_ILLUM_ARRAYS
  Type *energy;
  Vector3 *stream;
  Matrix3 *impuls;
#endif

  device_host_ptr_t() : illum(nullptr), int_scattering(nullptr), divstream(nullptr), divimpuls(nullptr),
                        normals(nullptr), areas(nullptr), volume(nullptr), directions(nullptr)
#ifdef ON_FULL_ILLUM_ARRAYS
                        ,
                        energy(nullptr), stream(nullptr), impuls(nullptr) {
  }
#endif
};

} // namespace geo
} // namespace cuda

extern cuda::geo::grid_directions_device_t *grid_dir_device;
extern cuda::geo::grid_device_t *grid_device;
extern cuda::geo::device_host_ptr_t device_host_ptr; ///< связь хоста с массивами внутри видеокарты (через эту структуры идёт обращение к массивам на карте)

#endif // CUDA_STRUCT_H