#ifdef USE_VTK
#include "global_types.h"
#include "utils.h"

#include "convert_vtk_geo.h"
#include "reader_vtk.h"
#include "writer_bin.h"

///\note: Формат хранения номеров граничных ячеек: i * CELL_SIZE + j

int FUNC_NAME(GetSurface)(int argc, char *argv[]) {

  if (argc != 3) {
    printf("Error input data!\n");
    printf("Input: ");
    printf("path\\file_vtk\n");
    printf("file_result.bin\\ \n");
    return e_completion_fail;
  }

  std::string file_vtk = argv[1];
  std::string file_result = argv[2];

  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  if (files_sys::vtk::Read(file_vtk, grid)) {
    return e_completion_fail;
  }

  std::vector<IntId> surface_id;
  GetBoundaryFacesId(grid, surface_id);

  std::vector<FaceCell> faces(surface_id.size());
  for (size_t i = 0; i < surface_id.size(); i++) {
    int id = surface_id[i];

    faces[i].face_id = id;
    for (vtkIdType k = 0; k < CELL_SIZE - 1; k++)
      faces[i].face[k] = Vector3(grid->GetCell(id / CELL_SIZE)->GetFace(id % CELL_SIZE)->GetPoints()->GetPoint(k));
  }

  return files_sys::bin::WriteSimple(file_result, faces);
}
#endif //! USE_VTK