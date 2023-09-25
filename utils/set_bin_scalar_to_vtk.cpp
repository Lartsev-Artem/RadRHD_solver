#include "global_def.h"
#include "prj_config.h"

#include "reader_bin.h"
#include "reader_txt.h"
#include "reader_vtk.h"
#include "writer_vtk.h"

#include <vtkDoubleArray.h>

typedef const std::string &file;
typedef int data_type;

#ifdef USE_VTK
static int WriteDataToGrid(file name_file_grid, file name_file_data, file name_file_output, file name_data) {
  vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  if (files_sys::vtk::Read(name_file_grid, unstructured_grid))
    RETURN_ERR("Error reading the file vtk\n");

  const size_t n = unstructured_grid->GetNumberOfCells();

  vtkSmartPointer<vtkDoubleArray> DataArray = vtkSmartPointer<vtkDoubleArray>::New();

  std::vector<data_type> vector_data;

  std::string formats = name_file_data.substr(name_file_data.find(".") + 1, 3);
  if (formats == "txt") {
    if (files_sys::txt::ReadSimple(name_file_data, vector_data))
      RETURN_ERR("");
  } else if (formats == "bin") {
    if (files_sys::bin::ReadSimple(name_file_data, vector_data))
      RETURN_ERR("");
  } else {
    RETURN_ERR("Error formats data file. use .bin or .txt\n");
  }

  if (vector_data.size() != n) {
    RETURN_ERR("Error size grid and data\n");
  }

  for (size_t i = 0; i < n; i++) {
    DataArray->InsertNextTuple1(vector_data[i]);
  }

  DataArray->SetName(name_data.c_str());
  unstructured_grid->GetCellData()->AddArray(DataArray);

  return files_sys::vtk::WriteVtkGrid(name_file_output, unstructured_grid);
}

int SetScalarDataVtkFromFile(int argc, char *argv[]) {
  if (argc < 4) {
    printf("Error input data!\n");
    printf("Input Format: path\\file.vtk,  path\\data.bin,  path\\outfile.vtk\n");
    printf("Additional parameters: name_field_data");
    printf("Data format: size a b c ... \n");
    return 1;
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

  return 0;
}
#endif // USE_VTK