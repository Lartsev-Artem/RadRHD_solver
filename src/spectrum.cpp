#include "mpi_ext.h"
#include "reader_json.h"
#include "solvers_struct.h"
#include "string_utils.h"

#include "illum_main.h"

int main(int argc, char *argv[]) {

#ifdef SPECTRUM
  MPI_START(argc, argv);

  int N = 1;
  std::string file_config = "/home/artem/projects/solver/config/directories_cfg.json";
  if (argc > 1)
    file_config = argv[1];

  if (files_sys::json::ReadStartSettings(file_config, glb_files, &_solve_mode, &_hllc_cfg))
    return e_completion_fail;

  // if (is_number(argv[1])) {
  //   N = std::stoi(argv[1]);
  // }

  // if (argc > 2 && is_number(argv[2])) {
  //   N = std::stoi(argv[2]);
  // }

  illum::RunSpectrumModule(N);

  MPI_END;
#else
  WRITE_LOG_ERR("the spectrum solver is not included in the build. Turn on define SPECTRUM in solvers_config.h\n");
  return 1;
#endif
  return 0;
}