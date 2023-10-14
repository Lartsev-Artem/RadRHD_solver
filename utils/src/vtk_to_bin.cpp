#ifdef USE_VTK
#include "global_types.h"
#include "utils.h"

#include "convert_vtk_geo.h"
#include "reader_vtk.h"

#include "writer_bin.h"

int FUNC_NAME(WriteBinFromVtk)(int argc, char *argv[]) {

  if (argc != 3) {
    printf("Error input data!\n");
    printf("Input: ");
    printf("path\\grid.vtk\n");
    printf("output_address\\ \n");
    return e_completion_fail;
  }

  std::string name_file_vtk = argv[1];
  std::string base_address = argv[2];

  vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  if (files_sys::vtk::Read(name_file_vtk, unstructured_grid))
    RETURN_ERR("Error reading the file vtk\n");

  const int size_grid = unstructured_grid->GetNumberOfCells();

  if (unstructured_grid->GetCellData()->GetNumberOfArrays() == 0)
    RETURN_ERR("grid hasn't data\n");

  std::vector<Type> s_data;
  std::vector<Vector3> v_data;
  std::vector<Matrix3> t_data;

  for (size_t i = 0; i < unstructured_grid->GetCellData()->GetNumberOfArrays(); i++) {
    std::string name_data(unstructured_grid->GetCellData()->GetArrayName(i));
    vtkDataArray *data = unstructured_grid->GetCellData()->GetScalars(name_data.c_str());

    int components = data->GetNumberOfComponents();
    int size = (data->GetSize() - components) / components;
    switch (components) {
    case 1: {
      s_data.resize(size);
      for (int i = 0; i < size; i++) {
        s_data[i] = data->GetTuple1(i);
      }

      files_sys::bin::WriteSimple(base_address + name_data + ".bin", s_data);
      break;
    }
    case 3: {
      v_data.resize(size);
      for (int i = 0; i < size; i++) {
        Vector3 buf;
        for (int j = 0; j < components; j++) {
          buf[j] = data->GetTuple3(i)[j];
        }

        v_data[i] = buf;
      }

      files_sys::bin::WriteSimple(base_address + name_data + ".bin", v_data);
      break;
    }
    case 9: {
      t_data.resize(size);
      for (int i = 0; i < size; i++) {
        Matrix3 buf;
        for (int j = 0; j < components; j++) {
          buf.data()[j] = data->GetTuple9(i)[j];
        }
        t_data[i] = buf;
      }

      files_sys::bin::WriteSimple(base_address + name_data + ".bin", t_data);
      break;
    }

    default:
      break;
    }
  }

  return e_completion_success;
}

#endif //! USE_VTK