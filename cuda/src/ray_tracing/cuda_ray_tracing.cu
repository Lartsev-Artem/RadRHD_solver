#ifdef USE_CUDA

#include "cuda_ray_struct.h"

__device__ static Type RayIntersectsTriangle(const Ray &ray, const Face &triangle, Vector3 &intersection) {

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

__global__ void RayTracing(const int M, const Ray *rays, const int N, const Face *triangles, Intersection *) {

  const int i = blockIdx.x * blockDim.x + threadIdx.x; //ячейка
  const int k = blockIdx.y * blockDim.y + threadIdx.y; //направление

  if (i >= N || k >= M) {
    return;
  }

  Intersection[k * N + i].dist = RayIntersectsTriangle(rays[k], triangles[i], Intersection[k * N + i].p);
}

#include <algorithm>
#include <solvers_struct.h>

int IllumSum(geo_face_t f, geo_cell_t *c) {

  //определить ячейку
  Type sign = f.n.dot(dir);

  if (sign > 0 && c->sign > 0) //одного знака наверное
  {
    get_cell[f.r]
  } else {
    get_cell[f.l]; // r/l может и наоборот
  }
}

int run() {
  std::vector<Ray> rays;
  std::vector<Face> faces;
  std::vector<Intersection> inters(cell * dir);

  // pragma parallel for
  int dir_i;
  int cnt[dir] = {0}; //число ячеек реально пересекавшие лучи
  auto{
    [dir_i](const Intersection &a, const Intersection &b) {
      cnt[dir_i] += (a.dist > 0);
      return a.dist < b.dist;
    }
  }

  for (dir_i = 0; dir_i < dir; dir_i++) {
    std::sort(inters.begin() + dir_i * N, inters.begin() + (dir_i + 1) * N, cmp);
  }

  std::vector<geo_cells_t> cells;
  double sum = 0;
  for (size_t i = 0; i < dir; i++) {
    for (size_t k = 0; k < cnt[i]; k++) {
      dir += IllumSum(faces, cells.data());
    }
  }

  return 0;
}

#endif //! USE_CUDA