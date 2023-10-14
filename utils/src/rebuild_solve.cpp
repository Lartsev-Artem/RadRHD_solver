#ifdef USE_VTK
#include "global_types.h"
#include "utils.h"

#include "reader_vtk.h"
#include "writer_vtk.h"

#include "set_vtk_data.h"

int FUNC_NAME(RebuildSolve)(int argc, char *argv[]) {
  if (argc != 4) {
    printf("Error input data!\n");
    printf("Input:\n");
    printf("name_file_vtk\n");
    printf("address_solve (like \"path\\file\")\n");
    printf("max_number_of_iter\n");
    return e_completion_fail;
  }

  const std::string name_file_vtk = argv[1];
  const std::string address_solve = argv[2];
  const int max_number_of_iter = std::stoi(argv[3]);

  for (int i = 0; i < max_number_of_iter; i++) {
    {
      vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
      if (files_sys::vtk::Read(name_file_vtk, grid)) {
        return e_completion_fail;
      }
      if (SetSolutionFromFileToVtk(address_solve + std::to_string(i), grid) == e_completion_success) {
        if (files_sys::vtk::WriteVtkGrid(address_solve + std::to_string(i) + ".vtk", grid, true) == e_completion_success) {
          printf("add grid %d\n", i);
        }
      }
    }
  }
  return e_completion_success;
}
#endif