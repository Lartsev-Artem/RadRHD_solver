#include <mpi.h>

int main(int argc, char* argv[]) {
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  printf("Rank: %d/%d\n", rank, size);

  MPI_Finalize();
  return 0;
}

#ifdef OST_BIN_READER
//#ifdef SOLVERS
int ReadGeometryGrid(const std::string& file_cells,
                     const std::string& file_faces, grid_t& grid) {
  ReadGridGeo(file_faces, grid.faces);
  ReadGridGeo(file_cells, grid.cells);

  // READ_FILE(file_faces.c_str(), grid.faces, geo);
  // READ_FILE(file_cells.c_str(), grid.cells, geo);
  grid.size = grid.cells.size();
  return 0;
}

//#if defined HLLC || defined RHLLC
template <typename T>
int ReadValueGrid(const T main_dir, grid_t& grid) {
  READ_FILE(std::string(main_dir + "phys_val.bin").c_str(), grid.cells,
            phys_val);
  return 0;
}
template <typename T>
int ReadHllcInit(const T file_init_value, std::vector<elem_t>& cells) {
  READ_FILE(std::string(file_init_value).c_str(), cells, phys_val);
  return 0;
}
#endif