#ifdef USE_CUDA
#include "cuda_ray_geo_intersect.h"
#include "cuda_ray_interface.h"
#include "ray_tracing_const.h"

namespace ray = cuda::ray_tracing;

static __device__ void sort(const int size, Type *array) {

  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size - 1; j++) {
      if (array[j] > array[j + 1]) {
        Type buf = array[j];     // создали дополнительную переменную
        array[j] = array[j + 1]; // меняем местами
        array[j + 1] = buf;      // значения элементов
      }
    }
  }
}

__device__ Type ray::Rosh(const Ray &ray, Type t) {

  const Vector3 a = ray.dir;
  const Vector3 S = ray.orig;

  // wolfram
  constexpr Type max_potential = -1.82547; // потенциал роша
  constexpr Type m = 0.8795180722891566;

  // /// увеличена масса донора с 0.1М до 0.22 для полного затмения полярных областей расчётной области
  // constexpr Type max_potential = -1.9242410534319847; // потенциал роша
  // constexpr Type m = 0.7684210526315789;

  // математические расчёты
  Type x1 = S[0] + a[0] * t;
  Type x2 = (S[1] + a[1] * t) * (S[1] + a[1] * t);
  Type x3 = (S[2] + a[2] * t) * (S[2] + a[2] * t);

  return -max_potential - ((x1 - m) * (x1 - m) + x2 + x3) * 0.5 -
         m / sqrt((x1 - 1) * (x1 - 1) + x2 + x3) -
         (1 - m) / sqrt(x1 * x1 + x2 + x3);
}

__device__ int ray::FindRoshRoots(const int n, const Type a, const Type b, const Ray &ray, Type *roots) {
  /* поиск пересечений кривой с осью ОХ */

  Type step = (b - a) / (n - 1);
  int size = 0;

  Type fi_1 = Rosh(ray, a + 0 * step);
  Type fi_2;
  for (int i = 1; i < n; ++i) {

    fi_2 = Rosh(ray, a + (i)*step);
    //если функция по разные стороны от ось ОХ => пересекла
    if (fi_1 * fi_2 < 0) {
      roots[size] = ((2 * a + (2. * i - 1) * step) / 2.);
      size++;
    }
    fi_1 = fi_2;
  }

  sort(size, roots);
  return size;
}

__device__ int ray::GetIntersectionWithRosh(const Ray &ray, Type *dist) {

  Type roots_rosh_equation[6];
  int number_of_roots = FindRoshRoots(1001, -3, 3, ray, roots_rosh_equation); // 1001 --- число узлов, [-3,3] --- интервал поиска
                                                                              //нет пересечения
  if (number_of_roots == 6 || number_of_roots == 4)                           // 4,6 --- связь с особенностью поверхности
  {

    Vector3 intersection_point = ray.orig + ray.dir * roots_rosh_equation[1];
    if (intersection_point[0] < 0.3) // 0.3 --- связь с положением точки зрения и поверхности (пересечение с поверхностью донера, а не аккретора)
    {
      *dist = roots_rosh_equation[1];
      return e_ray_intersect_rosh;
    }
  }

  *dist = -1;
  return e_ray_intersect_none;
}

__device__ void ray::IntersectionWithPlane(const Face &face, const Ray &ray, Vector3 &result) {

  Type a, b, c, d; // параметры уравнения плоскости
  Type t;

  a = face.A[1] * (face.B[2] - face.C[2]) + face.B[1] * (face.C[2] - face.A[2]) + face.C[1] * (face.A[2] - face.B[2]);
  b = face.A[0] * (face.C[2] - face.B[2]) + face.B[0] * (face.A[2] - face.C[2]) + face.C[0] * (face.B[2] - face.A[2]);
  c = face.A[0] * (face.B[1] - face.C[1]) + face.B[0] * (face.C[1] - face.A[1]) + face.C[0] * (face.A[1] - face.B[1]);
  d = face.A[0] * (face.C[1] * face.B[2] - face.B[1] * face.C[2]) + face.B[0] * (face.A[1] * face.C[2] - face.C[1] * face.A[2]) + face.C[0] * (face.B[1] * face.A[2] - face.A[1] * face.B[2]);

  t = -(a * ray.orig[0] + b * ray.orig[1] + c * ray.orig[2] + d) / (a * ray.dir[0] + b * ray.dir[1] + c * ray.dir[2]);

  result = ray.dir * t + ray.orig; // точка пересечения луча   с плоскостью!!!
}

__device__ static const Type c_x = ray_tracing::k_center_x;
__device__ static const Type c_y = ray_tracing::k_center_y;
__device__ static const Type c_z = ray_tracing::k_center_z;
__device__ static const Type radius_sphere = ray_tracing::k_radius_sphere;
__device__ static const Type R_in_2 = ray_tracing::k_internal_radius_disk * ray_tracing::k_internal_radius_disk;
__device__ static const Type R_ex_2 = ray_tracing::k_external_radius_disk * ray_tracing::k_external_radius_disk;

__device__ int ray::GetIntersectionWithSphereOfDisk(const Ray &ray) {

  Face plane_disk; // точки задающие плоскость диска (Wolfram)
  plane_disk.A = Vector3(1, 0, 0);
  plane_disk.B = Vector3(0, 0.9928768384869221, 0.11914522061843064);
  plane_disk.C = Vector3(2, 0, 0);

  Vector3 center(c_x, c_y, c_z);

  // базисные вектора, задающие наклон аккреционного диска вне расчетной плоскости (Wolfram)
  Vector3 vec_1(1, 0, 0);
  Vector3 vec_2(0, -0.992877, -0.119145);

  Vector3 res(0, 0, 0); // пересечние луча с плоскостью диска
  IntersectionWithPlane(plane_disk, ray, res);

  Vector3 local_intersection(0, 0, 0); // точка пересечения в локальных координатах плоскости
  for (int k = 0; k < 3; k++) {        // в плоскости
    local_intersection[0] += (res[k]) * vec_1[k];
    local_intersection[1] += (res[k]) * vec_2[k];
  }
  Type loc_x = local_intersection[0] - center[0];
  Type loc_y = local_intersection[1] - center[1];
  Type loc_square_of_distance = loc_x * loc_x + loc_y * loc_y;

  const Vector3 &n = ray.dir;
  Type A = n.dot(n);
  Type proj = (ray.orig - center).dot(n);

  // подкоренное выражение из уравнения пересечения сферы и прямой
  Type radical = 4. * (proj * proj - A * (1 - 2 * ray.orig[0] + ray.orig.dot(ray.orig) - radius_sphere * radius_sphere));

  if (radical >= 0) // есть пересечения со сферой
  {

    Type t = (n[0] - n.dot(ray.orig) - 0.5 * sqrt(radical)) / A;
    Vector3 in_sphere = n * t + ray.orig; // точка на сфере

    // не пересекает плоскость диска (в локальных координатах (без наклона))
    if (loc_square_of_distance <= R_in_2 || loc_square_of_distance >= R_ex_2) {
      return e_ray_intersect_sphere;
    }

    // с чем луч встречается раньше?
    Type distance_to_spehere = (ray.orig - in_sphere).norm();
    Type distance_to_plane = (ray.orig - res).norm();

    if (distance_to_spehere > distance_to_plane) {
      return e_ray_intersect_disk;
    } else {
      return e_ray_intersect_sphere;
    }
  } else //нет пересечения со сферой но внутри диска
  {
    if (loc_square_of_distance < R_ex_2 && (loc_square_of_distance > R_in_2)) {
      return e_ray_intersect_disk;
    }
  }

  return e_ray_intersect_none;
}

#endif