#ifdef MAKE_TRACE
#include "trace_main.h"

#include "trace_nodes.h"

#include "global_types.h"
#include "global_value.h"
#include "intersections.h"
#include "mpi_ext.h"
#include "reader_bin.h"
#include "reader_txt.h"
#include "solvers_struct.h"
#include "trace_to_face.h"
#include "writer_bin.h"

#include <chrono>
namespace tick = std::chrono;

int trace::RunTracesModule() {

  WRITE_LOG("Start RunTracesModule\n");

  int np = get_mpi_np();
  int myid = get_mpi_id();

  //-----------файлы с данными сетки. Построены здесь на метки USE_VTK-----------------------
  const std::string name_file_normals = glb_files.base_address + F_NORMALS;
  const std::string name_file_cells = glb_files.base_address + F_TRACE_GRID;
  const std::string name_file_vertex = glb_files.base_address + F_TRACE_VERTEX;
  const std::string name_file_neigh = glb_files.base_address + F_NEIGHBOR;

  //--------------------------------создающиеся файлы----------------------------------------

#ifdef TRANSFER_CELL_TO_FACE
  const std::string name_file_graph_bound = glb_files.graph_address + F_GRAPH_BOUND_FACE;
  const std::string name_file_graph_body = glb_files.graph_address + F_GRAPH_BODY_FACE;
#else
  const std::string name_file_state_face = glb_files.illum_geo_address + F_STATE_FACE;
#endif

  const std::string name_file_x = glb_files.illum_geo_address + F_X;
  const std::string name_file_x0_loc = glb_files.illum_geo_address + F_X0_LOC;

  // make
  std::vector<Face> grid;
  std::vector<Matrix4> vertexs;
  std::vector<Normals> normals;
  std::vector<IntId> neighbours;

  grid_directions_t grid_direction;

  auto start_clock = tick::steady_clock::now();

  uint32_t err = 0;
  err |= files_sys::bin::ReadSimple(glb_files.name_file_neigh, neighbours);
  err |= files_sys::bin::ReadSimple(name_file_cells, grid);
  err |= files_sys::bin::ReadSimple(name_file_vertex, vertexs);
  err |= files_sys::bin::ReadNormals(name_file_normals, normals);
  err |= files_sys::txt::ReadSphereDirectionCartesian(glb_files.name_file_sphere_direction, grid_direction);

#ifdef TRANSFER_CELL_TO_FACE
  grid_t geo_grid;
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_cells, geo_grid.cells);
  err |= files_sys::bin::ReadGridGeo(glb_files.name_file_geometry_faces, geo_grid.faces);
#endif

  WRITE_LOG("Reading time trace prebuild %lf\n", (double)tick::duration_cast<tick::milliseconds>(tick::steady_clock::now() - start_clock).count() / 1000.);

  const int count_cells = vertexs.size();
  std::vector<BasePointTetra> vec_x(count_cells);
  if (GetInterpolationNodes(vertexs, vec_x)) {
    D_LD;
  }

  if (myid == 0) {
    if (files_sys::bin::WriteSimple(name_file_x, vec_x)) {
      RETURN_ERR("Error writing %s\n", name_file_x.c_str());
    }
  }

  std::vector<IntId> sorted_graph;

  int num_cell;
  Vector3 direction;

  // массивы записи в файл:
  std::vector<cell_local> vec_x0;

#ifdef TRANSFER_CELL_TO_FACE
  std::vector<IntId> graph_bound_faces;
  std::vector<graph_pair_t> graph_cell_faces;
#else
  std::vector<bits_flag_t> face_states(count_cells, 0); //битовое поле: 0=> выходящая грань,  1=> входящая
#endif

  /*---------------------------------- далее FOR по направлениям----------------------------------*/
#if defined ONLY_ONE_DIRECTION
  for (int num_direction = 0; num_direction < 1; num_direction++)
#else
  for (int num_direction = myid; num_direction < grid_direction.size; num_direction += np)
#endif
  {

#if defined ONLY_ONE_DIRECTION
    num_direction = 0;
#endif

    direction = grid_direction.directions[num_direction].dir;

    DIE_IF(files_sys::bin::ReadSimple(glb_files.graph_address + F_GRAPH + std::to_string(num_direction) + ".bin", sorted_graph));

    vec_x0.clear();
    vec_x0.reserve(CELL_SIZE * count_cells);
#ifdef TRANSFER_CELL_TO_FACE
    graph_bound_faces.clear();
    graph_bound_faces.reserve(10000);
    graph_cell_faces.clear();
    graph_cell_faces.reserve(2 * geo_grid.cells.size());
#endif

    /*---------------------------------- далее FOR по ячейкам----------------------------------*/
    for (int h = 0; h < count_cells; ++h) {
      num_cell = sorted_graph[h];

      bits_flag_t face_state = 0;
      intersection::FindInAndOutFaces(direction, normals[num_cell], face_state);
#ifndef TRANSFER_CELL_TO_FACE
      face_states[num_cell] = face_state;
#endif

      for (ShortId num_out_face = 0; num_out_face < CELL_SIZE; ++num_out_face) {
        if (CHECK_BIT(face_state, num_out_face) == e_face_type_out) // выходящие грани
        {
          GetLocNodes(num_cell, num_out_face, grid, vertexs[num_cell],
                      face_state, direction, normals[num_cell], neighbours,
                      vec_x[num_cell], vec_x0);

#ifdef TRANSFER_CELL_TO_FACE
          graph_pair_t buf;
          buf.cell = num_cell;
          buf.loc_face = num_out_face;
          graph_cell_faces.push_back(buf); //это вместо graph
#endif
        }
#ifdef TRANSFER_CELL_TO_FACE
        //==e_face_type_in
        else if (neighbours[num_cell * CELL_SIZE + num_out_face] < 0) ///< сосед к текущей грани
        {
          graph_bound_faces.push_back(geo_grid.cells[num_cell].geo.id_faces[num_out_face]); //номера граней на сетке
        }
#endif
      }
    }
    /*---------------------------------- конец FOR по ячейкам----------------------------------*/

#ifdef TRANSFER_CELL_TO_FACE

    if (files_sys::bin::WriteSimple(name_file_graph_bound + std::to_string(num_direction) + ".bin", graph_bound_faces))
      RETURN_ERR("Error graph_bound_faces");

    if (files_sys::bin::WriteSimple(name_file_graph_body + std::to_string(num_direction) + ".bin", graph_cell_faces))
      RETURN_ERR("Error graph_cell_faces");

#ifndef DEBUG
    std::remove((glb_files.graph_address + F_GRAPH + std::to_string(num_direction) + ".bin").c_str()); //удаляем граф по ячейкам
#endif

#else
    if (files_sys::bin::WriteSimple(name_file_state_face + std::to_string(num_direction) + ".bin", face_states))
      RETURN_ERR("Error face_states");
#endif
#ifdef TRANSFER_CELL_TO_FACE
    WRITE_FILE_ELEM((name_file_x0_loc + "_s" + std::to_string(num_direction) + ".bin").c_str(), vec_x0, s);
    std::vector<face_loc_id_t> compress_id(vec_x0.size() / NODE_SIZE);
    face_loc_id_t val;
    for (size_t i = 0; i < compress_id.size(); i++) {
      val.a = vec_x0[i * NODE_SIZE + 0].in_face_id;
      val.b = vec_x0[i * NODE_SIZE + 1].in_face_id;
      val.c = vec_x0[i * NODE_SIZE + 2].in_face_id;
      compress_id[i].bits = val.bits;
    }

    files_sys::bin::WriteSimple(name_file_x0_loc + "_id" + std::to_string(num_direction) + ".bin", compress_id);
#else
    if (files_sys::bin::WriteSimple(name_file_x0_loc + std::to_string(num_direction) + ".bin", vec_x0))
      RETURN_ERR("Error vec_x0");
#endif

    WRITE_LOG("End trace direction number # %d\n", num_direction);
  }
  /*---------------------------------- конец FOR по направлениям----------------------------------*/

  WRITE_LOG("Full trace time: %lf\n", (double)tick::duration_cast<tick::milliseconds>(tick::steady_clock::now() - start_clock).count() / 1000.);

  MPI_BARRIER(MPI_COMM_WORLD);

  WRITE_LOG("End RunTracesModule\n");
  return e_completion_success;
}

#endif ///! MAKE_TRACE