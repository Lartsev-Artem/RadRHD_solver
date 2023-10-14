#include "ray_tracing_build_plane.h"
#include "ray_tracing_const.h"

#include "intersections.h"

#ifdef USE_VTK
#include "global_value.h"
#include "reader_bin.h"
#include "writer_vtk.h"

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkQuad.h>

void ray_tracing::MakeVtkPlane(vtkSmartPointer<vtkUnstructuredGrid> &image_plane) {

  // компоненты картинной плоскости
  vtkSmartPointer<vtkPoints> points_quad = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> quad_array = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

  const Vector3 angle_of_plane(-(k_width_plane / 2), -(k_height_plane / 2), 0); // угол плоскости. От него начинается заполнение всей плоскости
  constexpr Type step_x = k_width_plane / k_pixels_width;                       // ширина пикселя
  constexpr Type step_y = k_height_plane / k_pixels_height;                     // высота пикселя

  size_t quad_number = 0;
  for (int i = 0; i < k_pixels_width; ++i) {
    for (int j = 0; j < k_pixels_height; ++j) {

      Vector3 orig_2d(angle_of_plane(0) + i * step_x, angle_of_plane(1) + j * step_y, 0); ///< центр нового пикселя на плоскости

      points_quad->InsertNextPoint(orig_2d[0] - step_x / 2., orig_2d[1] - step_y / 2., 0);
      points_quad->InsertNextPoint(orig_2d[0] + step_x / 2., orig_2d[1] - step_y / 2., 0);
      points_quad->InsertNextPoint(orig_2d[0] + step_x / 2., orig_2d[1] + step_y / 2., 0);
      points_quad->InsertNextPoint(orig_2d[0] - step_x / 2., orig_2d[1] + step_y / 2., 0);

      quad->GetPointIds()->SetId(0, quad_number++);
      quad->GetPointIds()->SetId(1, quad_number++);
      quad->GetPointIds()->SetId(2, quad_number++);
      quad->GetPointIds()->SetId(3, quad_number++); // (№вершины, №точки)
      quad_array->InsertNextCell(quad);
    }
  }

  image_plane->SetPoints(points_quad);         // массив точек
  image_plane->SetCells(VTK_QUAD, quad_array); // массив ячеек

  return;
}

int ray_tracing::BuildVtkFromBin(const int number_of_planes, const std::string &files_plane) {

  vtkSmartPointer<vtkUnstructuredGrid> grid_plane = vtkSmartPointer<vtkUnstructuredGrid>::New();
  MakeVtkPlane(grid_plane);

  vtkSmartPointer<vtkDoubleArray> Illum_array = vtkSmartPointer<vtkDoubleArray>::New();
  std::vector<Type> energy_plane;

  for (int i = 0; i < number_of_planes; i++) {
    if (files_sys::bin::ReadSimple(files_plane + std::to_string(i) + ".bin", energy_plane))
      return e_completion_fail;

    D_L;
    Illum_array->SetArray(energy_plane.data(), energy_plane.size(), 1);
    D_L;
    grid_plane->GetCellData()->SetScalars(Illum_array);
    D_L;
    if (files_sys::vtk::WriteVtkGrid(files_plane + std::to_string(i) + ".vtk", grid_plane))
      return e_completion_fail;
  }

  return e_completion_success;
}

#endif

void ray_tracing::MakeRays(int num_frame, std::vector<Ray_t> &rays) {

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
      Vector3 orig_3d = basis.transpose() * orig_2d + center_ray.orig;                    // переход к 3d

      rays[i * k_pixels_height + j] = Ray_t(orig_3d, center_ray.direction);
    }
  }

  return;
}
