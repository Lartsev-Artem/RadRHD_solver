#ifdef USE_CUDA
#include "ray_tracing_main.h"
#include "cuda_ray_interface.h"
#include "global_value.h"
#include "ray_tracing_build_plane.h"
#include "ray_tracing_calc_illum.h"
#include "ray_tracing_const.h"

#include "reader_bin.h"
#include "writer_bin.h"
#include "writer_vtk.h"

int ray_tracing::RunRayTracing(const std::string &file_energy) {

  if (get_mpi_id() != 0) {
    return e_completion_success;
  }
#ifdef USE_VTK
  {
    vtkSmartPointer<vtkUnstructuredGrid> grid;
    MakeVtkPlane(grid);
    files_sys::vtk::WriteVtkGrid(glb_files.trace_address + F_IMAGE_PLANE + ".vtk", grid);
  }
#endif

  std::vector<FaceCell> faces;
  if (files_sys::bin::ReadSimple(glb_files.base_address + F_SURFACE, faces))
    return e_completion_fail;

  std::vector<Ray_t> rays(k_pixels_width * k_pixels_height);
  std::vector<IntId> intersections(k_pixels_width * k_pixels_height, -1);

  cuda::ray_tracing::interface::InitDevice(faces, rays);

  for (int i = 0; i < k_number_of_frame; i++) {
    MakeRays(i, rays);
    cuda::ray_tracing::interface::StartTracing(rays, intersections);
    files_sys::bin::WriteSimple(glb_files.trace_address + F_RAY_TRACE + std::to_string(i) + ".bin", intersections);
  }

  cuda::ray_tracing::interface::ClearDevice();

  WRITE_LOG("Ray Tracing ready\n");

  MakeEnergyAndCurve(file_energy);

  WRITE_LOG("Image plane ready\n");

  return e_completion_success;
}
#endif ///! USE_CUDA