#ifndef READ_NETGEN_H
#define READ_NETGEN_H
#include "utils.h"

#include "geo_types.h"
#include "global_def.h"

namespace utils {

template <typename Point, typename Cell>
int ReadNetgenMeshGrid(const std::string &name_file_in, std::vector<Point> &point, std::vector<Cell> &cell) {

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

template <typename Point, typename Cell>
int ReadNetgenGrid(const std::string &name_file_in, std::vector<Point> &point, std::vector<Cell> &cell, const size_t size = 3) {

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

} // namespace utils
#endif //! READ_NETGEN_H