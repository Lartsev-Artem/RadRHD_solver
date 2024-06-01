#include "mpi_ext.h"
#include "radRHD_main.h"
#include "reader_json.h"
#include "solvers_struct.h"

int main(int argc, char *argv[]) {

  // MPI_START(argc, argv);
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
#if defined RAD_RHD && defined SEPARATE_GPU && !defined SPECTRUM

  std::string file_config = "/home/artem/projects/solver/config/directories_cfg.json";
  if (argc > 1)
    file_config = argv[1];

  if (files_sys::json::ReadStartSettings(file_config, glb_files, &_solve_mode, &_hllc_cfg))
    return e_completion_fail;

  rad_rhd::RunRadRHDMpiModule();

#else
  WRITE_LOG_ERR("the rhllc solver is not included in the build. Use define RHLLC and SOLVER,SEPARATE_GPU\n");
#endif
  MPI_END;
  return 0;
}