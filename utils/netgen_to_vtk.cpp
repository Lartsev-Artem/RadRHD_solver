#include "global_def.h"
#include <string>
#include <vector>

namespace utils {

static int ReadNetgenGrid(const std::string &name_file_in, std::vector<Eigen::Vector3d> &point, std::vector<Eigen::Vector4i> &cell, const size_t size = 3) {

  std::ifstream ifile;
  OPEN_FSTREAM(ifile, name_file_in.c_str());

  size_t n;
  ifile >> n;
  point.resize(n);

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < size; j++) {
      ifile >> point[i][j];
    }
  }

  ifile >> n;
  cell.resize(n);
  for (size_t i = 0; i < n; i++) {
    ifile >> cell[i][0];
    for (size_t j = 0; j < size + 1; j++) {
      ifile >> cell[i][j];
    }
  }
  ifile.close();

  return e_completion_success;
}

static int WriteVtkFile(const std::string &name_file_out, const std::vector<Vector3> &point, const std::vector<Eigen::Vector4i> &cell, const size_t size = 3) {
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

int ReBuildNetgenToVTK(int argc, char *argv[]) {
  std::string name_file_in = "";
  std::string name_file_out = "";
  int size = 0;

  if (argc <= 3) {
    RETURN_ERR("Error input data!\n Input: path\\input_file, path\\output_file.vtk\n, size_task[2 or 3]");
  } else {
    name_file_in = argv[1];
    name_file_out = argv[2];
    size = std::stoi(argv[3]);
  }

  if (size != 2 && size != 3) {
    RETURN_ERR("bad size task. use 2d or 3d\n");
  }

  std::vector<Vector3> point;
  std::vector<Eigen::Vector4i> cell;

  ReadNetgenGrid(name_file_in, point, cell, size);
  WriteVtkFile(name_file_out, point, cell, size);

  return e_completion_success;
}
} // namespace utils