#include "reader_txt.h"

int ReadSphereDirectionСartesian(const std::string& file_sphere_direction,
                                 grid_directions_t& grid_direction) {
  std::ifstream ifile;
  OPEN_FSTREAM(ifile, file_sphere_direction.c_str());

  int N = 0;
  ifile >> N;

  grid_direction.size = N;
  grid_direction.directions.resize(N);

  for (int i = 0; i < N; i++) {
    ifile >> grid_direction.directions[i].area;
    ifile >> grid_direction.directions[i].dir[0] >>
        grid_direction.directions[i].dir[1] >>
        grid_direction.directions[i].dir[2];
  }
  ifile >> grid_direction.full_area;
  ifile.close();

  return e_completion_success;
}
