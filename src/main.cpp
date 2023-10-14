#include "graph_main.h"
#include "illum_main.h"
#include "mpi_ext.h"
#include "reader_json.h"
#include "trace_main.h"

#include <unistd.h>

#include "solvers_struct.h"

#include "ray_tracing_main.h"
#include "rhllc_main.h"

int main(int argc, char *argv[]) {

  MPI_START(argc, argv);

  files_sys::json::ReadStartSettings("/home/artem/projects/solver/config/directories_cfg.json",
                                     glb_files, &_solve_mode, &_hllc_cfg);

#ifdef ILLUM
  // graph::RunGraphModule();
  // trace::RunTracesModule();

  // if (get_mpi_id() == 0) {
  //   illum::RunIllumModule();
  // }
// GDB_ATTACH;
#endif

  ray_tracing::RunRayTracing(glb_files.solve_address + "0" + F_ENERGY);

  //  rhllc::RunRhllcModule();

  MPI_END;
  return 0;
}