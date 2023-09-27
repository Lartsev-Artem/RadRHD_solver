#ifdef USE_VTK
#include "write_data_to_vtk.h"

int FUNC_NAME(SetScalarDataVtkFromBinFile)(int argc, char *argv[]) {
  if (argc < 4) {
    printf("Error input data!\n");
    printf("Input Format: path\\file.vtk,  path\\data.bin,  path\\outfile.vtk\n");
    printf("Additional parameters: name_field_data");
    printf("Data format: size a b c ... \n");
    return e_completion_fail;
  }

  std::string file_grid = argv[1];
  std::string file_data = argv[2];
  std::string file_out = argv[3];
  std::string name_data = "data";

  if (argc > 4) {
    name_data = argv[4];
  }

  if (WriteDataToGrid(file_grid, file_data, file_out, name_data)) {
    RETURN_ERR("Error write data to vtk grid\n");
  }

  return e_completion_success;
}
#endif // USE_VTK