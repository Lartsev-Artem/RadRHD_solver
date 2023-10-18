#include "illum_calc_gpu_async.h"
#include "illum_init_data.h"
#include "ray_tracing_build_plane.h"
#include "ray_tracing_const.h"

#include "illum_calc_cpu.h"

#include "cuda_ray_interface.h"

#include "reader_bin.h"
#include "writer_bin.h"

using namespace illum;
using namespace ray_tracing;

#define F_DIRECTION_INTERPOLATION ""

int GetDirectionReInterpolation(int add_directions, const grid_directions_t &grid_dir) {

  std::vector<Face> direction_surface;
  if (files_sys::bin::ReadSimple("", direction_surface))
    return e_completion_fail;

  DIE_IF(grid_dir.size != direction_surface.size()); //разные сферы направлений

  std::vector<FaceCell> faces(direction_surface.size());
  for (size_t i = 0; i < direction_surface.size(); i++) {
    faces[i].face = direction_surface[i];
    faces[i].face_id = i;
  }
  direction_surface.clear();

  std::vector<Ray_t> rays(add_directions);
  std::vector<int> dir_id(add_directions);

  for (size_t i = 0; i < add_directions; i++) {

    Type start_angle = 2 * PI / add_directions * i; // текущий поворот картинной плоскости
    // выполняется поворот относительно центра масс на угол start_angle относительно нулевого положения
    const Vector3 start_ray = Vector3(cos(start_angle), sin(start_angle), k_height_above_center) + k_center_of_mass;
    const Vector3 end_ray(k_center_of_mass);

    rays[i].direction = -(end_ray - start_ray).normalized();
    rays[i].orig = Vector3::Zero(); //сфера направлений задана в нуле
  }

  cuda::ray_tracing::interface::InitDevice(faces, rays);
  cuda::ray_tracing::interface::TracingGrid(rays, dir_id);
  cuda::ray_tracing::interface::ClearDevice();

  grid_directions_t grid_plane;
  for (size_t i = 0; i < add_directions; i++) {
    grid_plane.directions[i].dir = rays[i].direction;
    grid_plane.directions[i].area = 1;
  }
  grid_plane.full_area = add_directions;
  grid_plane.size = add_directions;
  // write sphere directions

  return files_sys::bin::WriteSimple(glb_files.base_address + F_DIRECTION_INTERPOLATION, dir_id); //номера ячеек переинтерполяции для сетки направлений
}

int Func(const grid_directions_t &grid_dir, grid_t &grid) {
  std::vector<int> dir_node;
  if (files_sys::bin::WriteSimple(glb_files.base_address + F_DIRECTION_INTERPOLATION, dir_node))
    return e_completion_fail;

  // graph_main(sphere_direction)
  // trace_main(sphere_direction)

  // main_illum, calc...
  // save->scattering_dir

  // read

  // illum::cpu::CalculateAdditionalIllum()

  return;
}