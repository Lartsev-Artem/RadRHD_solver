
#include "ray_tracing_calc_illum.h"

#include "global_value.h"
#include "ray_tracing_const.h"

#include "reader_bin.h"
#include "writer_bin.h"
#include "writer_txt.h"

#ifdef USE_VTK
#include "reader_vtk.h"
#include "writer_vtk.h"
#include <vtkDoubleArray.h>

int ray_tracing::MakeEnergyAndCurve(const std::string &file_energy) {

  std::vector<IntId> intersections;
  std::vector<Type> energy;

  files_sys::bin::ReadSimple(file_energy, energy);

  vtkSmartPointer<vtkUnstructuredGrid> grid_plane = vtkSmartPointer<vtkUnstructuredGrid>::New();
  if (files_sys::vtk::Read(glb_files.trace_address + F_IMAGE_PLANE + ".vtk", grid_plane))
    return e_completion_fail;
  size_t size_plane = grid_plane->GetNumberOfCells();

  vtkSmartPointer<vtkDoubleArray> Illum_array = vtkSmartPointer<vtkDoubleArray>::New();
  Illum_array->SetNumberOfTuples(size_plane);

  std::vector<Type> curve_light(k_number_of_frame);

  for (int frame = 0; frame < k_number_of_frame; frame++) {

    curve_light[frame] = 0;

    if (files_sys::bin::ReadSimple(glb_files.trace_address + F_RAY_TRACE + std::to_string(frame) + ".bin", intersections))
      return e_completion_fail;

    for (size_t i = 0; i < intersections.size(); i++) {

      switch (intersections[i]) {
      case e_ray_intersect_none:
        Illum_array->SetTuple1(i, 0);
        break;

      case e_ray_intersect_rosh:
        Illum_array->SetTuple1(i, k_rosh_energy);
        break;
      case e_ray_intersect_disk:
        Illum_array->SetTuple1(i, k_disk_energy);
        break;

      case e_ray_intersect_sphere:
        Illum_array->SetTuple1(i, k_accretion_energy);
        break;

      default:
        DIE_IF(intersections[i] < 0);
        Illum_array->SetTuple1(i, energy[intersections[i] / CELL_SIZE]);
        break;
      }
      curve_light[frame] += Illum_array->GetTuple1(i);
    }

    grid_plane->GetCellData()->SetScalars(Illum_array);
    if (files_sys::vtk::WriteVtkGrid(glb_files.trace_address + F_IMAGE_PLANE + std::to_string(frame) + ".vtk", grid_plane))
      return e_completion_fail;
  }

  return files_sys::txt::WriteSimple(glb_files.trace_address + F_GLOSS_CURVE, curve_light);
}

#else

int ray_tracing::MakeEnergyAndCurve(const std::string &file_energy) {

  std::vector<IntId> intersections;
  std::vector<Type> energy;

  files_sys::bin::ReadSimple(file_energy, energy);

  std::vector<Type> Illum_array;
  std::vector<Type> curve_light(k_number_of_frame);

  for (int frame = 0; frame < k_number_of_frame; frame++) {

    curve_light[frame] = 0;
    if (files_sys::bin::ReadSimple(glb_files.trace_address + F_RAY_TRACE + std::to_string(frame) + ".bin", intersections))
      return e_completion_fail;

    Illum_array.assign(intersections.size(), 0);

    for (size_t i = 0; i < intersections.size(); i++) {

      switch (intersections[i]) {
      case e_ray_intersect_none:
        Illum_array[i] = 0;
        break;

      case e_ray_intersect_rosh:
        Illum_array[i] = k_rosh_energy;
        break;
      case e_ray_intersect_disk:
        Illum_array[i] = k_disk_energy;
        break;

      case e_ray_intersect_sphere:
        Illum_array[i] = k_accretion_energy;
        break;

      default:
        DIE_IF(intersections[i] < 0);
        Illum_array[i] = energy[intersections[i] / CELL_SIZE];
        break;
      }

      curve_light[frame] += Illum_array[i];
    }

    if (files_sys::bin::WriteSimple(glb_files.trace_address + F_IMAGE_PLANE + std::to_string(frame) + ".bin", Illum_array))
      return e_completion_fail;
  }

  return files_sys::txt::WriteSimple(glb_files.trace_address + F_GLOSS_CURVE, curve_light);
}

#endif //! USE_VTK
int ray_tracing::MakeIllumAndCurve(const std::string &base_file_illum) {

  std::vector<IntId> intersections;

  std::vector<Type> illum;
  std::vector<Type> Illum_array;
  std::vector<Type> curve_light(k_number_of_frame);

  for (int frame = 0; frame < k_number_of_frame; frame++) {

    files_sys::bin::ReadSimple(base_file_illum + std::to_string(frame) + ".bin", illum);

    curve_light[frame] = 0;
    if (files_sys::bin::ReadSimple(glb_files.trace_address + F_RAY_TRACE + std::to_string(frame) + ".bin", intersections))
      return e_completion_fail;

    Illum_array.assign(intersections.size(), 0);

    for (size_t i = 0; i < intersections.size(); i++) {

      switch (intersections[i]) {
      case e_ray_intersect_none:
        Illum_array[i] = 0;
        break;

      case e_ray_intersect_rosh:
        Illum_array[i] = k_rosh_energy;
        break;
      case e_ray_intersect_disk:
        Illum_array[i] = k_disk_energy;
        break;

      case e_ray_intersect_sphere:
        Illum_array[i] = k_accretion_energy;
        break;

      default:
        DIE_IF(intersections[i] < 0);
        Illum_array[i] = illum[intersections[i] / CELL_SIZE];
        break;
      }

      curve_light[frame] += (Illum_array[i] * (k_width_plane / k_pixels_width) * (k_height_plane / k_pixels_height));
    }

    if (files_sys::bin::WriteSimple(glb_files.trace_address + F_IMAGE_PLANE + std::to_string(frame) + ".bin", Illum_array))
      return e_completion_fail;
  }

  return files_sys::txt::WriteSimple(glb_files.trace_address + F_GLOSS_CURVE, curve_light);
}