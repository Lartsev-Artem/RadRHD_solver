#ifndef WRITE_METIS_H
#define WRITE_METIS_H
#include "utils.h"

#include "geo_types.h"
#include "global_def.h"

namespace utils {

template <typename Cell>
int WriteMetis(const std::string &name_file_out, const std::vector<Cell> &cell, const size_t size = 3) {
  std::ofstream ofile;
  OPEN_FSTREAM(ofile, name_file_out.c_str());

  ofile << cell.size() << "\n";
  for (size_t i = 0; i < cell.size(); i++) {
    for (size_t j = 0; j < size + 1; j++) {
      ofile << cell[i][j] << ' ';
    }
    ofile << '\n';
  }

  ofile.close();
  return e_completion_success;
}
} // namespace utils
#endif //! WRITE_METIS_H