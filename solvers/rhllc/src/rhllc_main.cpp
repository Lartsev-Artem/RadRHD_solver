
#if defined SOLVERS
#include "rhllc_main.h"

#include "reader_bin.h"
#include "writer_bin.h"

#include "solvers_struct.h"

int rhllc::RunRhllcModule() {
  grid_t grid;

  uint32_t err = 0;

  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, grid.faces);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, grid.cells);
  if (err) {
    RETURN_ERR("Error reading \n");
  }
  grid.InitMemory(grid.cells.size(), 0);

  /*
  run rhllc
  */

  files_sys::bin::WriteSolution(glb_files.solve_address + "0", grid);

  return e_completion_success;
}

#endif //! SOLVERS