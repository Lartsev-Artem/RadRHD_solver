#include "global_def.h"
#include "utils.h"

#include "geo_types.h"

#include <vector>

static int ReadNetgenGrid(const std::string &name_file_in, std::vector<Eigen::Vector3d> &point, std::vector<Eigen::Vector4i> &cell) {

  std::ifstream ifile;
  OPEN_FSTREAM(ifile, name_file_in.c_str());

  std::string buf = "";

  while (strcmp(buf.c_str(), "volumeelements")) {
    std::getline(ifile, buf);
  }
  int N, a;
  ifile >> N;

  cell.resize(N);
  for (size_t i = 0; i < N; i++) {
    ifile >> a >> a >> cell[i][0] >> cell[i][1] >> cell[i][2] >> cell[i][3];
  }

  while (strcmp(buf.c_str(), "points")) {
    std::getline(ifile, buf);
  }
  ifile >> N;

  point.resize(N);
  for (size_t i = 0; i < N; i++) {
    ifile >> point[i][0] >> point[i][1] >> point[i][2];
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

int FUNC_NAME(ReBuildNetgenMeshToVTK)(int argc, char *argv[]) {
  std::string name_file_in = "";
  std::string name_file_out = "";

  if (argc <= 2) {
    RETURN_ERR("Error input data!\n Input: path\\input_file, path\\output_file.vtk\n");
  } else {
    name_file_in = argv[1];
    name_file_out = argv[2];
  }

  std::vector<Vector3> point;
  std::vector<Eigen::Vector4i> cell;

  ReadNetgenGrid(name_file_in, point, cell);
  WriteVtkFile(name_file_out, point, cell, 3);

  return e_completion_success;
}
