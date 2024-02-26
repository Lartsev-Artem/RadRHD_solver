#include "global_def.h"
#include "utils.h"

#include "geo_types.h"
#include "read_netgen.h"
#include "write_metis.h"

#include <vector>

int FUNC_NAME(NetgenMeshToMetis)(int argc, char *argv[]) {
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

  utils::ReadNetgenMeshGrid(name_file_in, point, cell);
  utils::WriteMetis(name_file_out, cell, 3);

  return e_completion_success;
}