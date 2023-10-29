#ifdef USE_VTK
#include "global_types.h"
#include "utils.h"

#include "convert_vtk_geo.h"
#include "reader_vtk.h"

#include <cctype>

int FUNC_NAME(MakeSeries1dProjection)(int argc, char *argv[]) {

  if (argc != 5) {
    printf("Error input data!\n");
    printf("Input: ");
    printf("path\\Solve_i\n");
    printf("base_address\\ \n");
    printf("Axis: X,Y,Z, R\n");
    printf("number of projections\n");
    return e_completion_fail;
  }

  std::string base_file_vtk = argv[1];
  std::string base_address = argv[2];
  char axis = toupper(argv[3][0]);
  const int max_iter = std::stoi(argv[4]);

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

  for (int iter = 0; iter < max_iter; iter++) {

    std::string file_vtk = base_file_vtk + std::to_string(iter) + ".vtk";

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
      OPEN_FSTREAM(ofile, (base_address + name_data + "1d_" + axis + "_" + std::to_string(iter) + ".txt").c_str());
      if (components == 1) {

        for (int i = 0; i < size_grid; i++)
          ofile << ((axis_idx == 4) ? centers[i].norm() : centers[i][axis_idx]) << ' ' << data->GetTuple1(i) << "\n";

      } else {
        VectorX v(components);
        for (int i = 0; i < size_grid; i++) {
          double *tuple = data->GetTuple(i);
          for (int j = 0; j < components; j++) {
            v[j] = tuple[j];
          }

          ofile << ((axis_idx == 4) ? centers[i].norm() : centers[i][axis_idx]) << ' ' << v.norm() << "\n";
        }
      }
      ofile.close();
    }
  }

  return e_completion_success;
}
#endif