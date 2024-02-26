#ifndef WRITE_FILE_VTK_H
#define WRITE_FILE_VTK_H
#include "geo_types.h"
#include "global_def.h"
#include "utils.h"

namespace utils {

template <typename file, typename Point, typename Cell>
int WriteVtkFile(const file &name_file_out, const std::vector<Point> &point, const std::vector<Cell> &cell, const size_t size = 3) {
  std::ofstream ofile;
  OPEN_FSTREAM(ofile, name_file_out.c_str());

  ofile << "# vtk DataFile Version 3.0\n";
  ofile << "Elemensts Volumes\n";
  ofile << "ASCII\n";
  ofile << "DATASET UNSTRUCTURED_GRID\n";
  ofile << "POINTS " << point.size() << " double\n";

  for (size_t i = 0; i < point.size(); i++) {
    ofile << std::setprecision(16) << point[i][0] << ' ' << point[i][1] << " " << point[i][2] << '\n';
  }

  ofile << "CELLS " << cell.size() << " " << cell.size() * (size + 2) << '\n';

  for (size_t i = 0; i < cell.size(); i++) {
    ofile << size + 1 << " ";
    for (size_t k = 0; k < size + 1; k++) {
      ofile << cell[i][k] - 1 << ' ';
    }
    ofile << '\n';
    // ofile << size + 1 << " " << cell[i][0] - 1 << ' ' << cell[i][1] - 1 << ' ' << cell[i][2] - 1 << ' ' << cell[i][3] - 1 << '\n';
  }

  int vtkType = 10; // VTK_TETRA
  if (size == 2) {
    vtkType = 5; // VTK_TRIANGLE
  }

  ofile << "CELL_TYPES " << cell.size() << "\n";
  for (size_t i = 0; i < cell.size(); i++) {
    ofile << vtkType << '\n';
  }

  ofile.close();
  return e_completion_success;
}
} // namespace utils
#endif //! WRITE_FILE_VTK_H