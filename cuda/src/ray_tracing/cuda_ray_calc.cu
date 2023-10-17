#if defined USE_CUDA

#include "cuda_ray_calc.h"
#include "cuda_ray_geo_intersect.h"
#include "cuda_ray_interface.h"
#include "ray_tracing_const.h"

namespace ray = cuda::ray_tracing;

__device__ Type cuda::ray_tracing::RayIntersectsTriangle(const Ray &ray, const Face &triangle, Vector3 &intersection) {

  constexpr Type eps = 1e-10;
  Vector3 edge1 = triangle.B - triangle.A;
  Vector3 edge2 = triangle.C - triangle.A;
  Vector3 h = ray.dir.cross(edge2);

  Type det = edge1.dot(h);

  if (det > -eps && det < eps)
    return -1; // луч параллелен плоскости

  Type inv_det = 1.0 / det;
  Vector3 s = ray.orig - triangle.A;
  Type u = inv_det * s.dot(h);

  if (u < 0.0 || u > 1.0)
    return -1;

  Vector3 q = s.cross(edge1);
  Type v = inv_det * ray.dir.dot(q);

  if (v < 0.0 || u + v > 1.0)
    return -1;

  Type dist = inv_det * edge2.dot(q); // расстояние до плоскости

  if (dist > 0.0) // ray intersection
  {
    intersection = ray.orig + ray.dir * dist;
    return dist;
  }

  return -1; //Это означает, что линия пересекла треугольник против хода луча
}

__global__ void cuda::ray_tracing::RayTracing(const int M, const Ray *rays, const int N, const Face *triangles, int *intersections) {

  const int dir = blockIdx.x * blockDim.x + threadIdx.x; //ячейка

  if (dir >= M) {
    return;
  }

  int loc_id = e_ray_intersect_none;
  Type loc_dist = 1e10;

  Type rosh_dist;
  if (GetIntersectionWithRosh(rays[dir], &rosh_dist) != e_ray_intersect_none) {
    loc_dist = rosh_dist;
    loc_id = e_ray_intersect_rosh;
  }

  //по всей границе
  Vector3 p;
  for (int i = 0; i < N; i++) {

    Type d = RayIntersectsTriangle(rays[dir], triangles[i], p);

    //новое пересечение есть и оно ближе имеющегося
    if (d > 0 && d < loc_dist) {
      loc_id = triangles[i].id;
      loc_dist = d;
    }
  }

  /*
    сюда пересечения с внешней геометрией
  */
  intersections[dir] = loc_id;
}

__global__ void cuda::ray_tracing::InnerRayTracing(const int M, const Ray *rays, const int N, const Face *triangles, int *intersections) {

  const int dir = blockIdx.x * blockDim.x + threadIdx.x; //ячейка

  if (dir >= M) {
    return;
  }

  int loc_id = e_ray_intersect_none;
  Type loc_dist = 1e10;

  loc_id = GetIntersectionWithSphereOfDisk(rays[dir]);

  if (loc_id != e_ray_intersect_none) {
    intersections[dir] = loc_id;
    return;
  }

  //по всей границе
  Vector3 p;
  for (int i = 0; i < N; i++) {

    Type d = RayIntersectsTriangle(rays[dir], triangles[i], p);

    //новое пересечение есть и оно ближе имеющегося
    if (d > 0 && d < loc_dist) {
      loc_id = triangles[i].id;
      loc_dist = d;
    }
  }

  /*
    сюда пересечения с внешней геометрией
  */
  intersections[dir] = loc_id;
}

#endif //! USE_CUDA