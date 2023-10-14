#ifdef USE_VTK
#ifndef WRITE_DATA_TO_VTK_H
#define WRITE_DATA_TO_VTK_H
#include "utils.h"

#include "global_def.h"
#include "reader_bin.h"
#include "reader_txt.h"
#include "reader_vtk.h"
#include "writer_vtk.h"

#include <vtkDoubleArray.h>

namespace utils {
typedef int data_type;

template <typename file>
int WriteDataToGrid(file name_file_grid, file name_file_data, file name_file_output, file name_data) {

  vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  if (files_sys::vtk::Read(name_file_grid, unstructured_grid))
    RETURN_ERR("Error reading the file vtk\n");

  const size_t n = unstructured_grid->GetNumberOfCells();

  vtkSmartPointer<vtkDoubleArray> DataArray = vtkSmartPointer<vtkDoubleArray>::New();

  std::vector<data_type> vector_data;

  std::string formats = name_file_data.substr(name_file_data.find(".") + 1, 3);
  if (formats == "txt") {
    if (files_sys::txt::ReadSimple(name_file_data, vector_data))
      return e_completion_fail;
  } else if (formats == "bin") {
    if (files_sys::bin::ReadSimple(name_file_data, vector_data))
      return e_completion_fail;
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
} // namespace utils
#endif //! WRITE_DATA_TO_VTK_H
#endif //#ifdef USE_VTK