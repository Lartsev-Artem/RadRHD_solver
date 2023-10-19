#include "graph_main.h"
#include "illum_main.h"
#include "mpi_ext.h"
#include "reader_json.h"
#include "trace_main.h"

#include <unistd.h>

#include "solvers_struct.h"

#include "illum_add_main.h"
#include "illum_add_prebuild.h"
#include "ray_tracing_const.h"
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
  graph::RunGraphModule();
  trace::RunTracesModule();

  illum::additional_direction::PreBuildAddDirections(ray_tracing::k_number_of_frame);

  if (files_sys::json::ReadStartSettings(file_config, glb_files, &_solve_mode, &_hllc_cfg))
    return e_completion_fail;

  if (get_mpi_np() == 1) {
    illum::RunIllumModule();
  } else {
    illum::RunIllumMpiModule();
  }

#endif

  illum::additional_direction::RunModule();

  // ray_tracing::RunRayTracing(glb_files.solve_address + "0" + F_ENERGY);
  ray_tracing::RunRayTracing(glb_files.add_dir_address + F_ILLUM);

  //  rhllc::RunRhllcModule();

  MPI_END;
  return 0;
}