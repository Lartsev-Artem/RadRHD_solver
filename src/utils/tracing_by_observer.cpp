#include "mpi_ext.h"
#include "reader_json.h"
#include "solvers_struct.h"

int TracingObserver();

int main(int argc, char *argv[]) {

#ifdef TRANSFER_CELL_TO_FACE
  MPI_START(argc, argv);

  std::string file_config = "/home/artem/projects/solver/config/directories_cfg.json";
  if (argc > 1)
    file_config = argv[1];

  if (files_sys::json::ReadStartSettings(file_config, glb_files, &_solve_mode, &_hllc_cfg))
    return e_completion_fail;

  TracingObserver();

  MPI_END;
#else
  D_LD;
#endif
  return 0;
}