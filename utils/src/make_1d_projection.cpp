#ifdef USE_VTK
#include "global_types.h"
#include "utils.h"

#include "convert_vtk_geo.h"
#include "reader_vtk.h"

#include <cctype>

int FUNC_NAME(Make1dProjection)(int argc, char *argv[]) {

  if (argc != 4) {
    printf("Error input data!\n");
    printf("Input: ");
    printf("path\\grid.vtk\n");
    printf("base_address\\ \n");
    printf("Axis: X,Y,Z, R\n");
    return e_completion_fail;
  }

  std::string file_vtk = argv[1];
  std::string base_address = argv[2];
  char axis = toupper(argv[3][0]);

  int axis_idx;
  switch (axis) {
  case 'X':
    axis_idx = 0;
    break;

  case 'Y':
    axis_idx = 1;
    break;

  case 'Z':
    axis_idx = 2;
    break;

  case 'R':
    axis_idx = 4;
    break;

  default:
    RETURN_ERR("unknow axis %c\n", axis);
  }

  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  if (files_sys::vtk::Read(file_vtk, grid)) {
    return e_completion_fail;
  }

  if (grid->GetCellData()->GetNumberOfArrays() == 0)
    RETURN_ERR("grid hasn't data\n");

  std::vector<Vector3> centers;
  GetCentersOfCells3D(grid, centers);

  const int size_grid = grid->GetNumberOfCells();

  for (int i = 0; i < grid->GetCellData()->GetNumberOfArrays(); i++) {

    std::string name_data(grid->GetCellData()->GetArrayName(i));

    vtkDataArray *data = grid->GetCellData()->GetScalars(name_data.c_str());

    int components = data->GetNumberOfComponents();
    int size = (data->GetSize() - components) / components;

    if (size != size_grid) {
      printf("bad size data. grid= %d, data= %d\n", size_grid, size);
      continue;
    }

    std::ofstream ofile;
    OPEN_FSTREAM(ofile, (base_address + name_data + "1d_" + axis + ".txt").c_str());
    if (components == 1) {
      if (axis_idx != 4) {
        for (int i = 0; i < size_grid; i++)
          ofile << centers[i][axis_idx] << ' ' << data->GetTuple1(i) << "\n";
      } else {
        for (int i = 0; i < size_grid; i++)
          ofile << centers[i].norm() << ' ' << data->GetTuple1(i) << "\n";
      }

    } else {
      VectorX v(components);
      for (int i = 0; i < size_grid; i++) {
        double *tuple = data->GetTuple(i);
        for (int j = 0; j < components; j++) {
          v[j] = tuple[j];
        }

        if (axis_idx != 4)
          ofile << centers[i][axis_idx] << ' ' << v.norm() << "\n";
        else
          ofile << centers[i].norm() << ' ' << v.norm() << "\n";
      }
    }
    ofile.close();
  }

  return e_completion_success;
}
#endif