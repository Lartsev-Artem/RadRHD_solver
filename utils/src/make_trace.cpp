#ifdef USE_VTK
#include "global_types.h"
#include "utils.h"

#include "convert_vtk_geo.h"
#include "reader_vtk.h"
#include "writer_txt.h"

#include "intersections.h"

int FUNC_NAME(MakeTrace)(int argc, char *argv[]) {

  if (argc != 9) {
    printf("Error input data!\n");
    printf("Input: ");
    printf("path\\file_vtk\n");
    printf("file_result.txt\\ \n");
    printf("direction: X Y Z\n");
    printf("origin: X Y Z\n");
    return e_completion_fail;
  }

  std::string file_vtk = argv[1];
  std::string file_result = argv[2];

  Vector3 dir(std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]));
  Vector3 orig(std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]));

  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  if (files_sys::vtk::Read(file_vtk, grid)) {
    return e_completion_fail;
  }

  std::vector<IntId> intersections;
  intersection::GetIntersectionCellId(Ray_t(orig, dir), grid, intersections);

  return files_sys::txt::WriteSimple(file_result, intersections);
}
#endif //! USE_VTK