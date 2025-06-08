#include "graph_main.h"
#include "mpi_ext.h"
#include "reader_json.h"
#include "solvers_struct.h"

int main(int argc, char *argv[]) {

#ifdef BUILD_GRAPH
  MPI_START(argc, argv);

  std::string file_config = "/home/artem/projects/solver/config/directories_cfg.json";
  if (argc > 1)
    file_config = argv[1];

  if (files_sys::json::ReadStartSettings(file_config, glb_files, &_solve_mode, &_hllc_cfg))
    return e_completion_fail;

  graph::RunGraphModule();
  MPI_END;
#else
  WRITE_LOG_ERR("the graph builder is not included in the build. Use define BUILD_GRAPH\n");
  return 1;
#endif
  return 0;
}