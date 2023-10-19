#include "writer_txt.h"

int files_sys::txt::WriteSphereDirectionCartesian(const std::string &file_sphere_direction, const grid_directions_t &grid_direction) {

  std::ofstream ofile;
  OPEN_FSTREAM(ofile, file_sphere_direction.c_str());

  ofile << grid_direction.size << "\n";

  for (int i = 0; i < grid_direction.size; i++) {
    ofile << std::setprecision(16) << grid_direction.directions[i].area << "\n";
    ofile << std::setprecision(16)
          << grid_direction.directions[i].dir[0] << ' '
          << grid_direction.directions[i].dir[1] << ' '
          << grid_direction.directions[i].dir[2] << "\n";
  }
  ofile << std::setprecision(16) << grid_direction.full_area;
  ofile.close();

  return e_completion_success;
}
