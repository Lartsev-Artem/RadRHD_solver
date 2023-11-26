#include "illum_main.h"
#include "mpi_ext.h"
#include "reader_json.h"
#include "solvers_struct.h"

int main(int argc, char *argv[]) {

#if defined SEPARATE_GPU
  MPI_START(argc, argv);

  std::string file_config = "/home/artem/projects/solver/config/directories_cfg.json";
  if (argc > 1)
    file_config = argv[1];

  if (files_sys::json::ReadStartSettings(file_config, glb_files, &_solve_mode, &_hllc_cfg))
    return e_completion_fail;

  if (get_mpi_np() == 1) {
    D_LD;
  }
#ifdef USE_MPI
  else {
    illum::RunIllumMultiGpuModule();
  }
#endif

  MPI_END;
#else
  WRITE_LOG_ERR("the illum solver is not included in the build. Use define ILLUM and SOLVER\n");
  return 1;
#endif

  return 0;
}