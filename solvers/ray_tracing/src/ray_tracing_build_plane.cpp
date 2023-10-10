#include "ray_tracing_build_plane.h"
#include "ray_tracing_const.h"

#include "intersections.h"

void ray_tracing::MakePlane(int num_frame, std::vector<Ray_t> &rays) {

  Type start_angle = 2 * PI / k_number_of_frame * num_frame; // текущий поворот картинной плоскости

  // выполняется поворот относительно центра масс на угол start_angle относительно нулевого положения
  const Vector3 start_ray = Vector3(cos(start_angle), sin(start_angle), k_height_above_center) + k_center_of_mass;
  const Vector3 end_ray(k_center_of_mass);

  Ray_t center_ray(start_ray, (end_ray - start_ray).normalized()); ///< центр плоскости

  Matrix3 basis; ///< локальный базис картинной плоскости
  intersection::SetBasis(center_ray.direction, basis);

  const Vector3 angle_of_plane(-(k_width_plane / 2), -(k_height_plane / 2), 0); // угол плоскости. От него начинается заполнение всей плоскости
  constexpr Type step_x = k_width_plane / k_pixels_width;                       // ширина пикселя
  constexpr Type step_y = k_height_plane / k_pixels_height;                     // высота пикселя

  rays.resize(k_pixels_width * k_pixels_height);
  //формируем всю плоскость
  for (int i = 0; i < k_pixels_width; ++i) {
    for (int j = 0; j < k_pixels_height; ++j) {

      Vector3 orig_2d(angle_of_plane(0) + i * step_x, angle_of_plane(1) + j * step_y, 0); ///< центр нового пикселя на плоскости
      Vector3 orig_3d = basis * orig_2d + center_ray.orig;                                // переход к 3d

      rays[i * k_pixels_height + j] = Ray_t(orig_3d, center_ray.direction);
    }
  }

  return;
}
