#ifdef USE_VTK
#include "utils.h"
#include "write_data_to_vtk.h"

int FUNC_NAME(AddScalarDataVtkFromBinFile)(int argc, char *argv[]) {
  if (argc < 3) {
    printf("Error input data!\n");
    printf("Input Format: path\\file.vtk,  path\\data.bin");
    printf("Additional parameters: name_field_data");
    printf("Data format: size a b c ... \n");
    return 1;
  }

  std::string file_grid = argv[1];
  std::string file_data = argv[2];
  std::string name_data = "data";

  if (argc > 3) {
    name_data = argv[3];
  }

  if (utils::WriteDataToGrid(file_grid, file_data, file_grid, name_data)) {
    RETURN_ERR("Error write data to vtk grid\n");
  }

  return 0;
}
#endif // USE_VTK