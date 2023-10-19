#ifdef USE_VTK
#include "utils.h"
#include "write_data_to_vtk.h"

int FUNC_NAME(SetScalarDataVtkFromBinFile)(int argc, char *argv[]) {
  if (argc < 5) {
    printf("Error input data!\n");
    printf("Input Format: path\\file.vtk,  path\\data.bin,  path\\outfile.vtk\n");
    printf("Input data: int, double");
    printf("Additional parameters: name_field_data");
    printf("Data format: size a b c ... \n");
    return e_completion_fail;
  }

  std::string file_grid = argv[1];
  std::string file_data = argv[2];
  std::string file_out = argv[3];
  std::string name_data = "data";
  const char type_idx = argv[3][0];

  if (argc > 5) {
    name_data = argv[5];
  }

  switch (type_idx) {
  case 'i':
    if (utils::WriteDataToGrid(file_grid, file_data, file_out, name_data, int(0))) {
      RETURN_ERR("Error write data to vtk grid\n");
    }

  case 'd':
    if (utils::WriteDataToGrid(file_grid, file_data, file_out, name_data, double(0))) {
      RETURN_ERR("Error write data to vtk grid\n");
    }
    break;

  default:
    RETURN_ERR("unknown type data %s\n", argv[3]);
  }

  return e_completion_success;
}
#endif // USE_VTK