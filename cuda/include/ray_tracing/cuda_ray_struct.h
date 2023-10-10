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

typedef Eigen::Vector3d Vector3;
typedef Eigen::Matrix3d Matrix3;
typedef double Type;

namespace cuda::ray_tracing {

struct Ray {
  Vector3 orig;
  Vector3 dir;
};

struct Face {
  int id;
  Vector3 A;
  Vector3 B;
  Vector3 C;
};

struct Intersection {
  int id;
  Type dist;
  Vector3 point;
  __device__ Intersection(const int i = -1, const Type d = -1, const Vector3 p = Vector3::Zero()) : id(i), dist(d), point(p) {}
};

} // namespace cuda::ray_tracing

#endif //! CUDA_RAY_STRUCT_H