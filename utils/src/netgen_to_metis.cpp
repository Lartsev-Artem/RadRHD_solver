#include "global_def.h"
#include "utils.h"

#include "geo_types.h"
#include "read_netgen.h"
#include "write_metis.h"

#include <vector>

int FUNC_NAME(NetgenToMetis)(int argc, char *argv[]) {
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

  utils::ReadNetgenGrid(name_file_in, point, cell, size);
  utils::WriteMetis(name_file_out, cell, size);

  return e_completion_success;
}