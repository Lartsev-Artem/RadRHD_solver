#include "mpi_ext.h"
#include "reader_json.h"
#include "solvers_struct.h"

#include "spec_all.h"

int main(int argc, char *argv[]) {

#ifdef SPECTRUM
  MPI_START(argc, argv);

  std::string file_config = "/home/artem/projects/solver/config/directories_cfg.json";
  if (argc > 1)
    file_config = argv[1];

  if (files_sys::json::ReadStartSettings(file_config, glb_files, &_solve_mode, &_hllc_cfg))
    return e_completion_fail;

  illum::spec::RunIllumSpectrumModule();

  MPI_END;
#else
  WRITE_LOG_ERR("the spectrum solver is not included in the build. Turn on define SPECTRUM in solvers_config.h\n");
  return 1;
#endif
  return 0;
}