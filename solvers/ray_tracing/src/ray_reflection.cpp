#ifdef USE_CUDA
#include "ray_tracing_main.h"

#include "global_value.h"
#include "linear_alg.h"

#include "cuda_ray_interface.h"

#include "reader_bin.h"
#include "reader_txt.h"

#include "writer_bin.h"
#include "writer_txt.h"

/**
 * @brief Создать серию файлов для каждого из направлений содержащих номера направлений до отражения и номер соответствующей грани
 *
 * @param[in] grid_dir  сфера направлений
 * @param[in] grid_surface поверхность сетки
 * @param[in] short_intersection интерполяция отраженного направления
 * @return err_code
 */
static int MakeOrderedListReflections(const grid_directions_t &grid_dir,
                                      const std::vector<FaceCell> &grid_surface,
                                      const std::vector<std::vector<int>> &intersection) {

  grid_t grid;
  if (files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, grid.faces))
    return e_completion_fail;
  if (files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, grid.cells))
    return e_completion_fail;

  std::vector<IntId> sorted_id_bound_face; ///< номера граней определенной границы

  for (size_t dir = 0; dir < grid_dir.size; dir++) {

    std::vector<int> dir_before_reflection;
    dir_before_reflection.reserve(2 * grid_surface.size());

    //определение переинтерполяции
    for (size_t i = 0; i < intersection.size(); i++) //всего направлений
    {
      for (size_t j = 0; j < intersection[i].size(); j++) //по каждой ячейке границы
      {
        if (intersection[i][j] == dir) //переотражение в текущее направление
        {
          int id = grid_surface[j].face_id;
          int num_face = grid.cells[id / CELL_SIZE].geo.id_faces[id % CELL_SIZE];
          dir_before_reflection.push_back(num_face); //ячейка для которой выполнено равенство
          dir_before_reflection.push_back(i);        //исходное направление
        }
      }
    }

    if (files_sys::bin::ReadSimple(glb_files.graph_address + F_GRAPH_BOUND_FACE + std::to_string(dir) + ".bin", sorted_id_bound_face))
      return e_completion_fail;

    std::vector<IntId> num_bound_faces;
    std::vector<uint16_t> id_direction;
    num_bound_faces.reserve(sorted_id_bound_face.size());
    id_direction.reserve(sorted_id_bound_face.size());

    //сохранение минимального набора данных для каждого направления
    for (auto num_face : sorted_id_bound_face) {
      for (size_t i = 0; i < dir_before_reflection.size(); i += 2) {
        if (num_face == dir_before_reflection[i]) {
          const IdType neigh_id = grid.faces[num_face].geo.id_r;
          if (neigh_id == e_bound_lock) //сохраняем только непрозрачную стенку
          {
            num_bound_faces.push_back(num_face);                  // num_face
            id_direction.push_back(dir_before_reflection[i + 1]); // dir_before_reflection
          }
        }
      }
    }

    files_sys::bin::WriteSimple(glb_files.illum_geo_address + F_DIRECTION_REFLECTION + std::to_string(dir) + ".bin", num_bound_faces);
    files_sys::bin::WriteSimple(glb_files.illum_geo_address + F_FACES_REFLECTION + std::to_string(dir) + ".bin", id_direction);
  }

  return e_completion_success;
}

int ray_tracing::FindReflectionDirections() {

  const std::string base_adr_output = glb_files.illum_geo_address;

  std::vector<Face> direction_surface;
  if (files_sys::bin::ReadSimple(glb_files.base_address + F_SURFACE_SPHERE_DIRECTION, direction_surface))
    return e_completion_fail;

  grid_directions_t grid_dir;
  if (files_sys::txt::ReadSphereDirectionCartesian(glb_files.base_address + F_DIRECTION_GRID, grid_dir))
    return e_completion_fail;

  DIE_IF(grid_dir.size != direction_surface.size()); //разные сферы направлений

  std::vector<FaceCell> dir_faces(direction_surface.size());
  for (size_t i = 0; i < direction_surface.size(); i++) {
    dir_faces[i].face = direction_surface[i];
    dir_faces[i].face_id = i;
  }
  direction_surface.clear();

  std::vector<FaceCell> grid_surface;
  if (files_sys::bin::ReadSimple(glb_files.base_address + F_SURFACE, grid_surface))
    return e_completion_fail;

  std::vector<Normals> normals;
  if (files_sys::bin::ReadNormals(glb_files.base_address + F_NORMALS, normals))
    return e_completion_fail;

  std::vector<Ray_t> rays(grid_surface.size());
  cuda::ray_tracing::interface::InitDevice(dir_faces, rays);

  std::vector<std::vector<IntId>> intersections(grid_dir.size);
  for (size_t i = 0; i < grid_dir.size; i++) {
    intersections[i].resize(grid_surface.size(), -1);
  }

  for (size_t num_dir = 0; num_dir < grid_dir.size; num_dir++) {
    const Vector3 &dir = grid_dir.directions[num_dir].dir;

    // Init rays
    for (size_t i = 0; i < rays.size(); i++) {
      int cell = grid_surface[i].face_id / CELL_SIZE;
      int face = grid_surface[i].face_id % CELL_SIZE;
      rays[i].direction = reflection_dir(dir, normals[cell].n[face]);
      rays[i].orig = Vector3::Zero(); //сфера направлений задана в начале координат

      // WRITE_LOG("%d dir= {%lf,%lf,%lf}\n", cell, rays[i].direction[0], rays[i].direction[1], rays[i].direction[2]);
    }

    cuda::ray_tracing::interface::StartTracingGrid(rays, intersections[num_dir]);
  }
  cuda::ray_tracing::interface::ClearDevice();

  normals.clear();
  rays.clear();
  dir_faces.clear();

  DIE_IF(MakeOrderedListReflections(grid_dir, grid_surface, intersections) != e_completion_success)
  return e_completion_success;
}

#endif