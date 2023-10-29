#ifdef USE_VTK
#include "global_types.h"
#include "utils.h"

#include "convert_vtk_geo.h"
#include "reader_vtk.h"

#include "intersections.h"

int FUNC_NAME(MakeSeries1dRayProjection)(int argc, char *argv[]) {

  if (argc != 11) {
    printf("Error input data!\n");
    printf("Input: ");
    printf("path\\Solve_i\n");
    printf("address result\\ \n");
    printf("direction: X Y Z\n");
    printf("origin: X Y Z\n");
    printf("max iter\n");
    printf("Axis projection: X,Y,Z,R\n");
    return e_completion_fail;
  }

  std::string base_file_vtk = argv[1];
  std::string base_address = argv[2];

  Vector3 dir(std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]));
  Vector3 orig(std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]));

  const int max_iter = std::stoi(argv[9]);

  char axis = toupper(argv[10][0]);

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

  std::vector<IntId> cells_id;
  std::vector<Vector3> centers;
  {
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    if (files_sys::vtk::Read(base_file_vtk + "0.vtk", grid)) {
      return e_completion_fail;
    }
    intersection::GetIntersectionCellId(Ray_t(orig, dir), grid, cells_id);
    GetCentersOfCells3D(grid, centers);
  }
  std::cout << "Size trace ray: " << cells_id.size() << "\n";

  for (int iter = 0; iter < max_iter; iter++) {

    std::string file_vtk = base_file_vtk + std::to_string(iter) + ".vtk";

    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    if (files_sys::vtk::Read(file_vtk, grid)) {
      return e_completion_fail;
    }

    if (grid->GetCellData()->GetNumberOfArrays() == 0)
      RETURN_ERR("grid hasn't data\n");

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

        for (int id = 0; id < cells_id.size(); id++) {
          int i = cells_id[id];
          ofile << ((axis_idx == 4) ? centers[i].norm() : centers[i][axis_idx]) << ' ' << data->GetTuple1(i) << "\n";
        }

      } else {
        VectorX v(components);
        for (int id = 0; id < cells_id.size(); id++) {
          int i = cells_id[id];
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