#ifdef USE_VTK
#include "utils.h"
#include "write_data_to_vtk.h"

int FUNC_NAME(AddScalarDataVtkFromBinFile)(int argc, char *argv[]) {
  if (argc < 4) {
    printf("Error input data!\n");
    printf("Input Format: path\\file.vtk,  path\\data.bin");
    printf("Input data: int, double");
    printf("Additional parameters: name_field_data\n");
    printf("Data format: size a b c ... \n");
    return 1;
  }

  std::string file_grid = argv[1];
  std::string file_data = argv[2];
  std::string name_data = "data";
  const char type_idx = argv[3][0];

  if (argc > 4) {
    name_data = argv[4];
  }

  switch (type_idx) {
  case 'i':
    if (utils::WriteDataToGrid(file_grid, file_data, file_grid, name_data, int(0))) {
      RETURN_ERR("Error write data to vtk grid\n");
    }
    break;

  case 'd':
    if (utils::WriteDataToGrid(file_grid, file_data, file_grid, name_data, double(0))) {
      RETURN_ERR("Error write data to vtk grid\n");
    }
    break;

  default:
    RETURN_ERR("unknown type data %s\n", argv[3]);
  }

  return e_completion_success;
}
#endif // USE_VTK