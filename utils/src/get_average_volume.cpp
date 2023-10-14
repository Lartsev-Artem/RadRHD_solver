#ifdef USE_VTK
#include "global_types.h"
#include "utils.h"

#include "convert_vtk_geo.h"
#include "reader_vtk.h"

int FUNC_NAME(GetAverageSize3D)(int argc, char *argv[]) {

  if (argc != 2) {
    printf("Error input data!\n");
    printf("Input: ");
    printf("path\\grid.vtk\n");
    return e_completion_fail;
  }

  std::string file_vtk = argv[1];

  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  if (files_sys::vtk::Read(file_vtk, grid)) {
    return e_completion_fail;
  }

  std::vector<double> volumes;
  GetVolume3D(grid, volumes);

  std::vector<Normals> norms;
  std::vector<double> areas;
  GetNormalAndAreas3D(grid, norms, areas);

  double min = 1e10;
  double mid = 0;
  for (size_t i = 0; i < volumes.size(); i++) {
    double surface_area = 0;
    for (int j = 0; j < CELL_SIZE; j++) {
      surface_area += areas[i * CELL_SIZE + j]; // площадь поверхности
    }

    double r = (3. / 4.) * (volumes[i]) / surface_area;
    mid += r;
    min = std::min(min, r);
  }
  mid /= volumes.size();

  printf("average size= %.16lf\n", mid);
  printf("min size= %.16lf\n", min);

  return e_completion_success;
}

#endif