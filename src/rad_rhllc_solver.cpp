#include "mpi_ext.h"
#include "radRHD_main.h"
#include "reader_json.h"
#include "solvers_struct.h"

int main(int argc, char *argv[]) {

#ifdef RAD_RHD
  MPI_START(argc, argv);

  std::string file_config = "/home/artem/projects/solver/config/directories_cfg.json";
  if (argc > 1)
    file_config = argv[1];

  if (files_sys::json::ReadStartSettings(file_config, glb_files, &_solve_mode, &_hllc_cfg))
    return e_completion_fail;

  // rad_rhd::RadRHD_ConstRadStateTest();
  // rad_rhd::RunRadRHDModule();
  rad_rhd::RunRadRHDMpiModule();

  MPI_END;
#else
  WRITE_LOG_ERR("the rhllc solver is not included in the build. Use define RHLLC and SOLVER\n");
  return 1;
#endif
  return 0;
}