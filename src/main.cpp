//#include <mpi.h>

// #include "make_internal_format.h"
// #include "reader_bin.h"

// #include "writer_bin.h"

#include "graph_main.h"
#include "mpi_ext.h"
#include "reader_json.h"

int SetScalarDataVtkFromFile(int argc, char *argv[]);

int main(int argc, char *argv[]) {

  MPI_START(argc, argv);

  // files_sys::json::ReadStartSettings("/home/artem/projects/solver/config/directories_cfg.json", glb_files);

  // BuildDataFromVTK(glb_files);

  // SetScalarDataVtkFromFile(argc, argv);

  /*
  graph test
  DataArray->SetNumberOfTuples(n);
  for (size_t i = 0; i < n; i++) {
    DataArray->SetTuple1(vector_data[i], i);
  }
  */
  graph::RunGraphModule();

  // return 0;

  // std::vector<int> a;
  // files_sys::bin::ReadSimple("foo.bin", a);
  // files_sys::bin::WriteSimple("build/foo.bin", a);

  MPI_END;
  return 0;
}

#ifdef OST_BIN_READER
//#ifdef SOLVERS
int ReadGeometryGrid(const std::string &file_cells,
                     const std::string &file_faces, grid_t &grid) {
  ReadGridGeo(file_faces, grid.faces);
  ReadGridGeo(file_cells, grid.cells);

  // READ_FILE(file_faces.c_str(), grid.faces, geo);
  // READ_FILE(file_cells.c_str(), grid.cells, geo);
  grid.size = grid.cells.size();
  return 0;
}

//#if defined HLLC || defined RHLLC
template <typename T>
int ReadValueGrid(const T main_dir, grid_t &grid) {
  READ_FILE(std::string(main_dir + "phys_val.bin").c_str(), grid.cells,
            phys_val);
  return 0;
}
template <typename T>
int ReadHllcInit(const T file_init_value, std::vector<elem_t> &cells) {
  READ_FILE(std::string(file_init_value).c_str(), cells, phys_val);
  return 0;
}
#endif