#ifdef SOLVERS
#include "illum_main.h"

#include "global_types.h"
#include "illum_calc_cpu.h"

#include "reader_bin.h"
#include "reader_txt.h"
#include "writer_bin.h"

#include "cuda_interface.h"

int illum::RunIllumModule() {

  grid_t grid;
  grid_directions_t grid_direction;

  std::vector<IntId> neighbours;

  std::vector<BasePointTetra> vec_x;
  std::vector<std::vector<State>> face_states;
  std::vector<std::vector<cell_local>> vec_x0;
  std::vector<std::vector<IntId>> sorted_id_cell;
  std::vector<Type> vec_res_bound;

  uint32_t err = 0;
  err |= files_sys::bin::ReadSimple(glb_files.name_file_neigh, neighbours);
  err |= files_sys::txt::ReadSphereDirectionСartesian(glb_files.name_file_sphere_direction, grid_direction);
  err |= files_sys::bin::ReadRadiationTrace(grid_direction.size, glb_files,
                                            vec_x, face_states, vec_x0, sorted_id_cell, vec_res_bound);

  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, grid.faces);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, grid.cells);

  if (err) {
    RETURN_ERR("Error reading \n");
  }

  grid.InitMemory(grid.cells.size(), grid_direction.size);

#ifdef USE_CUDA
  cuda::interface::InitDevice(glb_files.base_address, grid_direction, grid, 0, grid_direction.size);
#endif

  cpu::CalculateIllum(grid_direction, face_states, neighbours,
                      vec_x0, vec_x, sorted_id_cell, grid);

  cpu::CalculateIllumParam(grid_direction, grid);

  files_sys::bin::WriteSolution(glb_files.solve_address + "0", grid);

#ifdef USE_CUDA
  cuda::interface::CudaWait();
  cuda::interface::ClearDevice();
#endif

  return e_completion_success;
}

#endif //! SOLVERS