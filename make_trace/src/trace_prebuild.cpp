#if defined USE_VTK
#include "trace_prebuild.h"
#include "global_value.h"

#include "convert_vtk_geo.h"
#include "geo_types.h"
#include "reader_vtk.h"
#include "writer_bin.h"

int trace::PreBuild(const global_files_t &glb_files) {
  //-----------файлы с данными сетки. Построены здесь на метки USE_VTK-----------------------
  const std::string name_file_cells = glb_files.base_address + F_TRACE_GRID;
  const std::string name_file_vertex = glb_files.base_address + F_TRACE_VERTEX;

  vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

  if (files_sys::vtk::Read(glb_files.name_file_vtk, unstructured_grid)) {
    RETURN_ERR("Error reading the file vtk\n");
  }

  std::vector<Face> faces;
  GetFacesPoints(unstructured_grid, faces);
  files_sys::bin::WriteSimple(name_file_cells, faces);
  faces.clear();

  const vtkIdType n = unstructured_grid->GetNumberOfCells();
  std::vector<Matrix4> vertexs(n);

  for (vtkIdType i = 0; i < n; i++) {
    GetVertexMatrix(i, unstructured_grid, vertexs[i]);
  }

  files_sys::bin::WriteSimple(name_file_vertex, vertexs);

  return e_completion_success;
}

#endif //! USE_VTK