/**
 * @file cuda_ray_struct.h
 * @brief Структуры модуля трассировки лучей
 *
 */

#if !defined CUDA_RAY_STRUCT_H && defined USE_CUDA
#define CUDA_RAY_STRUCT_H

// workaround issue between gcc >= 4.7 and cuda 5.5 (совместимость с компиляторами (см. док. Eigen3: https://eigen.tuxfamily.org/dox/TopicCUDA.html))
#if (defined __GNUC__) && (__GNUC__ > 4 || __GNUC_MINOR__ >= 7)
#undef _GLIBCXX_ATOMIC_BUILTINS
#undef _GLIBCXX_USE_INT128
#endif
#include <Eigen/Dense> /// требует компиляции с флагом --expt-relaxed-constexpr

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

/*! \addtogroup ray_tracing Модуль трассировки лучей, построения картинной плоскости и кривых блеска
    @{
*/

typedef Eigen::Vector3d Vector3;
typedef Eigen::Matrix3d Matrix3;
typedef double Type;

namespace cuda::ray_tracing {

/**
 * @brief луч
 *
 */
struct Ray {
  Vector3 orig; ///< начало луча
  Vector3 dir;  ///< направление
};

/**
 * @brief Грань. Соответствует структуре ::FaceCell
 *
 */
struct Face {
  int id;    ///<  номер грани в глобальной нумерации
  Vector3 A; ///< вершина треугольника
  Vector3 B; ///<вершина треугольника
  Vector3 C; ///<вершина треугольника
};

/**
 * @brief Структура пересечения
 *
 */
struct Intersection {
  int id;        ///<  номер ячейки пересечения
  Type dist;     ///< расстояние до пересечения
  Vector3 point; ///< точка пересечения
  __device__ Intersection(const int i = -1, const Type d = -1, const Vector3 p = Vector3::Zero()) : id(i), dist(d), point(p) {}
};

} // namespace cuda::ray_tracing

#endif //! CUDA_RAY_STRUCT_H