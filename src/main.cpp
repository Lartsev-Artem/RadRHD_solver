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

  std::string file_config = "/home/artem/projects/solver/config/directories_cfg.json";
  if (argc > 1)
    file_config = argv[1];

  if (files_sys::json::ReadStartSettings(file_config, glb_files, &_solve_mode, &_hllc_cfg))
    return e_completion_fail;

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