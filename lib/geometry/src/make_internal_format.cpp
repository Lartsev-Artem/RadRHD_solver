#include "make_internal_format.h"
#ifdef USE_VTK

#include "reader_bin.h"
#include "reader_vtk.h"
#include "writer_bin.h"
#include "writer_txt.h"

#include "convert_vtk_geo.h"

static int SetTypeOfBound(const std::vector<Vector3> &centers,
                          const std::vector<Normals> &normals, std::vector<int> &neighbors) {
  const size_t n = centers.size();
  for (size_t num_cell = 0; num_cell < n; num_cell++) {
    Vector3 P = centers[num_cell];
    for (size_t num_face = 0; num_face < CELL_SIZE; num_face++) {
      const size_t id = num_cell * CELL_SIZE + num_face;

      if (neighbors[id] < 0) // граница
      {
#if GEOMETRY_TYPE == Cone_JET
        if ((normals[num_cell].n[num_face] - Vector3(-1, 0, 0)).norm() < 1e-7) {
          if (Vector2(P[1], P[2]).norm() < 0.01)
            neighbors[id] = e_bound_out_source; // ИСТОЧНИК джет
          else
            neighbors[id] = e_bound_free; // свободная поверхность

        } else if ((normals[num_cell].n[num_face] - Vector3(1, 0, 0)).norm() < 1e-7) {
          neighbors[id] = e_bound_free; // свободная поверхность
        } else {
          neighbors[id] = e_bound_lock; // боковая поверхность
        }
#endif

#if GEOMETRY_TYPE == Cone
        if ((normals[num_cell].n[num_face] - Vector3(-1, 0, 0)).norm() < 1e-7) {
          neighbors[id] = e_bound_out_source; // излучающее дно
        } else if ((normals[num_cell].n[num_face] - Vector3(1, 0, 0)).norm() < 1e-7) {
          neighbors[id] = e_bound_free; // свободная поверхность
        } else {
          neighbors[id] = e_bound_lock; // боковая поверхность
        }
#endif

#if GEOMETRY_TYPE == Sphere
        if ((P - kCenterPoint).norm() > kInternalRadius)
          neighbors[id] = e_bound_free; // внешняя сфера
        else
          neighbors[id] = e_bound_inner_source; // внутренняя сфера
#endif

#if GEOMETRY_TYPE == Cube
        neighbors[id] = e_bound_free;
#endif

#if GEOMETRY_TYPE == Step
        if ((normals[num_cell].n[num_face] - Vector3(-1, 0, 0)).norm() < 1e-5) {
          neighbors[id] = e_bound_free; // левая стенка
        } else if ((normals[num_cell].n[num_face] - Vector3(1, 0, 0)).norm() < 1e-5 && P[0] > 0.7) {
          neighbors[id] = e_bound_free; // горизонтальные стенки
        } else {
          neighbors[id] = e_bound_lock; //ступенька
        }

#endif
#if GEOMETRY_TYPE == Cylinder
        if ((normals[num_cell].n[num_face] - Vector3(-1, 0, 0)).norm() < 1e-5 &&
            Vector2(P[1], P[2]).norm() < 0.2) {
          neighbors[id] = e_bound_out_source; // основание
        } else {
          neighbors[id] = e_bound_free; // боковая поверхность и др. основание
        }

#endif
      } // if
    }   // face
  }     // cell

  return e_completion_success;
}

static int ReWriteNeighByType(const std::string &name_file_pairs, const std::string &name_file_normals, const std::string &name_file_centers) {

  std::vector<int> neighbors;
  std::vector<Normals> normals;
  std::vector<Vector3> centers;

  files_sys::bin::ReadSimple(name_file_centers, centers);
  files_sys::bin::ReadSimple(name_file_pairs, neighbors);
  files_sys::bin::ReadNormals(name_file_normals, normals);

  SetTypeOfBound(centers, normals, neighbors);

  return files_sys::bin::WriteSimple(name_file_pairs, neighbors);
}

static int WriteNeighbor(const std::string &name_file_neighbors, const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid) {

  std::vector<int> neighbors;

#if NUMBER_OF_MEASUREMENTS == 3
  if (GetNeighborFace3D(unstructured_grid, neighbors))
#elif NUMBER_OF_MEASUREMENTS == 2
  if (GetNeighborFace2D(unstructured_grid, neighbors))
#endif
    RETURN_ERR("Neighbor wasn't found\n");

  return files_sys::bin::WriteSimple(name_file_neighbors, neighbors);
}

static int WriteNormalAndAreas(const std::string &name_file_normals, const std::string name_file_areas, vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid) {

  std::vector<Normals> normals;
  std::vector<Type> areas;
#if NUMBER_OF_MEASUREMENTS == 3
  if (GetNormalAndAreas3D(unstructured_grid, normals, areas))
#elif NUMBER_OF_MEASUREMENTS == 2
  if (GetNormalAndAreas3D(unstructured_grid, normals, areas))
#endif
    RETURN_ERR("normals wasn't found\n");

  if (files_sys::bin::WriteSimple(name_file_areas, areas) == e_completion_success) {
    return files_sys::bin::WriteNormals(name_file_normals, normals);
  }

  return e_completion_fail;
}

static int WriteVolume(const std::string name_file_volume, vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid) {

  std::vector<Type> volumes;

#if NUMBER_OF_MEASUREMENTS == 3
  if (GetVolume3D(unstructured_grid, volumes))
#elif NUMBER_OF_MEASUREMENTS == 2
  if (GetVolume2D(unstructured_grid, volumes))
#endif
    RETURN_ERR("volume wasn't found\n");

  return files_sys::bin::WriteSimple(name_file_volume, volumes);
}

static int WriteCentersOfCells(const std::string name_file_centers, const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid) {

  std::vector<Vector3> centers;

#if NUMBER_OF_MEASUREMENTS == 3
  if (GetCentersOfCells3D(unstructured_grid, centers))
#elif NUMBER_OF_MEASUREMENTS == 2
  if (GetCentersOfCells2D(unstructured_grid, centers))
#endif
    RETURN_ERR("centers cells wasn't found\n");

  return files_sys::bin::WriteSimple(name_file_centers, centers);
}

#if NUMBER_OF_MEASUREMENTS == 3

static int WriteCentersOfFaces(const std::string &name_file_centers_faces, const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid) {

  std::vector<Vector3> centers_faces;

  if (GetCentersOfFaces3D(unstructured_grid, centers_faces)) {
    RETURN_ERR("centers faces wasn't found\n");
  }

  return files_sys::bin::WriteSimple(name_file_centers_faces, centers_faces);
}

static int WriteInitBoundarySetCells(const std::string file_boundary_cells,
                                     const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid) {
  std::set<IntId> boundary_cells;

  if (GetBoundaryCells(unstructured_grid, boundary_cells)) {
    RETURN_ERR("boundary cells wasn't found\n");
  }

  return files_sys::txt::WriteSimple(file_boundary_cells, boundary_cells);
}
static int WriteInitBoundarySetInFaces(const std::string file_boundary_id_faces, const std::string file_boundary_faces,
                                       const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid) {

  std::set<IntId> inter_boundary_faces;
  if (GetInterBoundaryFacesOfSphere(unstructured_grid, inter_boundary_faces)) {
    RETURN_ERR("internal boundary faces wasn't found\n");
  }

  std::map<IntId, FaceCell> faces_coord;
  if (GetFacesById(inter_boundary_faces, unstructured_grid, faces_coord)) {
    RETURN_ERR("the faces are not rebuilt\n");
  }

  if (files_sys::txt::WriteSimple(file_boundary_id_faces, inter_boundary_faces)) {
    return e_completion_fail;
  }
  return files_sys::txt::WriteSimple(file_boundary_faces, faces_coord);
}

#endif // 3d

int BuildDataFromVTK(const global_files_t &glb_files) {

  vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  if (files_sys::vtk::Read(glb_files.name_file_vtk, unstructured_grid)) {
    RETURN_ERR("Error reading file vtk\n");
  }

  if (WriteNeighbor(glb_files.base_address + F_NEIGHBOR, unstructured_grid)) {
    RETURN_ERR("Error writing file %s\n", F_NEIGHBOR);
  }

  if (WriteNormalAndAreas(glb_files.base_address + F_NORMALS, glb_files.base_address + F_AREAS, unstructured_grid)) {
    RETURN_ERR("Error writing file %s\n", F_NORMALS);
  }

  if (WriteVolume(glb_files.base_address + F_VOLUME, unstructured_grid)) {
    RETURN_ERR("Error writing file %s\n", F_VOLUME);
  }

  if (WriteCentersOfCells(glb_files.base_address + F_CENTERS, unstructured_grid)) {
    RETURN_ERR("Error writing file %s\n", F_CENTERS);
  }

#if NUMBER_OF_MEASUREMENTS == 3

  if (WriteInitBoundarySetCells(glb_files.base_address + F_INIT_BOUND, unstructured_grid)) {
    RETURN_ERR("Error writing file %s\n", F_INIT_BOUND);
  }

  if (WriteInitBoundarySetInFaces(glb_files.base_address + F_INTERNAL_BOUND, glb_files.base_address + F_FACE_ID, unstructured_grid)) {
    RETURN_ERR("Error writing file %s or %s\n", F_INTERNAL_BOUND, F_FACE_ID);
  }

  if (WriteCentersOfFaces(glb_files.base_address + F_CENTER_FACE, unstructured_grid)) {
    RETURN_ERR("Error writing file %s\n", F_CENTER_FACE);
  }

#endif // 3d

  if (ReWriteNeighByType(glb_files.base_address + F_NEIGHBOR, glb_files.base_address + F_NORMALS,
                         glb_files.base_address + F_CENTERS)) {
    RETURN_ERR("Error rewriting file %s\n", F_NEIGHBOR);
  }

  return e_completion_success;
}

#endif //! USE_VTK