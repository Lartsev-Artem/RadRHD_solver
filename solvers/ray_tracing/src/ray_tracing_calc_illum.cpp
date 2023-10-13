#ifdef USE_VTK
#include "ray_tracing_calc_illum.h"

#include "global_value.h"
#include "ray_tracing_const.h"
#include "reader_bin.h"
#include "reader_vtk.h"
#include "writer_bin.h"
#include "writer_txt.h"
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

  for (int i = 0; i < k_number_of_frame; i++) {

    curve_light[i] = 0;

    if (files_sys::bin::ReadSimple(glb_files.trace_address + F_RAY_TRACE + std::to_string(i) + ".bin", intersections))
      return e_completion_fail;

    for (size_t i = 0; i < intersections.size(); i++) {

      switch (intersections[i]) {
      case e_ray_intersect_none:
        Illum_array->SetTuple1(i, 0);
        break;

      case e_ray_intersect_rosh:
        Illum_array->SetTuple1(i, 1);
        break;
      case e_ray_intersect_disk:
        Illum_array->SetTuple1(i, 5);
        break;

      case e_ray_intersect_sphere:
        Illum_array->SetTuple1(i, 10);
        break;

      default:
        DIE_IF(intersections[i] < 0);
        Illum_array->SetTuple1(i, energy[intersections[i]]);
        break;
      }
    }

    curve_light[i] += Illum_array->GetTuple1(i);
    grid_plane->GetCellData()->SetScalars(Illum_array);
    if (files_sys::vtk::WriteVtkGrid(glb_files.trace_address + F_IMAGE_PLANE + std::to_string(i) + ".vtk", grid_plane))
      return e_completion_fail;
  }

  return files_sys::txt::WriteSimple(glb_files.trace_address + F_GLOSS_CURVE, curve_light);
}

#endif //! USE_VTK