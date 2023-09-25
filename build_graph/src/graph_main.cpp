#if defined BUILD_GRAPH
#include "graph_main.h"

#include "dbgdef.h"
#include "graph_config.h"
#include "mpi_ext.h"
#include "reader_bin.h"
#include "reader_txt.h"
#include "writer_bin.h"

#include "graph_calc.h"
#include "graph_init_state.h"
#include "graph_struct.h"

#include <omp.h>

namespace graph {
std::vector<boundary_trace_t> bound_trace; ///< данные перетрассировки луча сквозь внутреннюю область
}
int graph::RunGraphModule() {

  WRITE_LOG("OMP num_threads: %d\n", omp_get_num_threads());

  int np = get_mpi_np();
  int myid = get_mpi_id();
  double t = -omp_get_wtime();

  std::vector<IntId> neighbours;          ///< соседние ячейки
  std::set<IntId> inter_boundary_face_id; ///< id внутренних граней [i * CELL_SIZE + j]
  std::vector<Normals> normals;           ///< нормали
  std::map<IntId, FaceCell> inter_faces;  ///< внутренние грани с ключом-номером ячейки
  grid_directions_t grid_dir;             ///< сфера направлений

  uint32_t err = 0;
  err |= files_sys::bin::ReadSimple(glb_files.name_file_neigh, neighbours);
  err |= files_sys::txt::ReadSimple(glb_files.base_address + F_INTERNAL_BOUND, inter_boundary_face_id);
  err |= files_sys::bin::ReadNormals(glb_files.base_address + F_NORMALS, normals);
  err |= files_sys::txt::ReadInitBoundarySetInFaces(glb_files.base_address + F_FACE_ID, inter_faces);
  err |= files_sys::txt::ReadSphereDirectionСartesian(glb_files.name_file_sphere_direction, grid_dir);

  if (err != 0) {
    RETURN_ERR("error during reading\n");
  }
  t += omp_get_wtime();

  WRITE_LOG("Inner boundary has %d faces\n", (int)inter_boundary_face_id.size());
  WRITE_LOG("Time reading in main proccess: %lf\n", t);

  t = -omp_get_wtime();

  const size_t num_cells = normals.size();
  WRITE_LOG("num_cells= %d\n", num_cells);

  std::vector<State> count_in_face(num_cells, 0);  ///< число входящих граней ячейки
  std::vector<State> count_def_face(num_cells, 0); ///< число определённых граней ячейки
  std::vector<IntId> graph(num_cells, 0);          ///< упорядоченный набор ячеек

  std::set<IntId> inner_part; ///< часть граничных ячеек через которые луч покидай расчётную область
  std::set<IntId> outer_part; ///< часть граничных ячеек через которые луч возвращается в расчётную область

  std::vector<State> faces_state; ///< состояние граней (определена не определена)

  bool flag = true;
  int count = 0;

#ifdef USE_OMP
  std::list<IntId> cur_el_OMP;         ///<текущая границы
  std::vector<IntId> next_step_el_OMP; ///< кандидаты на следующую границу
#else
  std::vector<IntId> cur_el;    ///<текущая границы
  std::set<IntId> next_step_el; ///< кандидаты на следующую границу
#endif // USE_OMP

#if defined ONLY_ONE_DIRECTION
  for (int cur_direction = 0; cur_direction < 1; cur_direction++)
#else
  for (int cur_direction = myid; cur_direction < grid_dir.size; cur_direction += np)
#endif // ONLY_ONE_DIRECTION
  {

    WRITE_LOG("Direction #: %d\n", cur_direction);

    flag = true;
    InitFacesState(neighbours, inter_faces, faces_state);
    Vector3 direction = grid_dir.directions[cur_direction].dir;

    DivideInnerBoundary(direction, normals, inter_boundary_face_id, inner_part, outer_part);

    bound_trace.clear();
    bound_trace.reserve(outer_part.size() * 3); // потенциально в массив могут войдут все ячейки внутренней границы

    //-------------------------------------
#ifdef USE_OMP
    FindNumberOfAllInnerFaceAndKnew(direction, normals, faces_state, count_in_face, count_def_face, next_step_el_OMP);
#else
    FindNumberOfAllInAndDefFaces(direction, normals, faces_state, count_in_face, count_def_face, next_step_el);
#endif // USE_OMP
       //-------------------------------------

    int count_graph = 0;         // число ячеек вошедших в граф
    graph.assign(num_cells, -1); // для отлавливания ошибочных направлений
    bool try_restart = true;

#ifdef USE_OMP
    while (next_step_el_OMP.size() && flag)
#else
    while (next_step_el.size())
#endif // USE_OMP
    {

#ifdef USE_OMP

#ifdef GRID_WITH_INNER_BOUNDARY
      IntId id_cell = FindCurCellWithHole(next_step_el_OMP, count_in_face, count_def_face, cur_el_OMP, inner_part, outter_part,
                                          inter_faces, neighbours, direction, normals);
#else
      IntId id_cell = FindCurCell(next_step_el_OMP, count_in_face, count_def_face, cur_el_OMP);
#endif // GRID_WITH_INNER_BOUNDARY

#else // no use omp

#ifdef GRID_WITH_INNER_BOUNDARY
      int cur_ret = FindCurFrontWithHole(direction, normals, inter_faces, next_step_el,
                                         count_in_face, inner_part, cur_el, count_def_face, outer_part);
#else
      int cur_ret = FindCurFront(next_step_el, count_in_face, count_def_face, cur_el);
#endif // GRID_WITH_INNER_BOUNDARY
#endif // USE_OMP

      if (cur_ret != e_completion_success) {
        WRITE_LOG("Warning proc: %d, dir= %d, processed %d cells", myid, cur_direction, count_graph);

        if (TryRestart(count_in_face, count_def_face, outer_part, cur_el, next_step_el) == e_completion_success) {
          WRITE_LOG("\n\n Warning!!! try_restart %d \n\n", cur_direction);
          continue;
        }
      }

#ifdef USE_OMP

      NewStep(neighbours, count_in_face, count_def_face, cur_el_OMP, next_step_el_OMP);

      for (auto el : cur_el_OMP) {
        graph[count_graph] = el;
        count_graph++;
      }
#else
      NewStep(neighbours, count_in_face, cur_el, count_def_face, next_step_el);

      for (auto el : cur_el) {
        graph[count_graph] = el;
        count_graph++;
      }
#endif // USE_OMP

    } // while

    DIE_IF_ACTION(count_graph < graph.size(), WRITE_LOG_ERR("Error size graph[%d] %d\n", cur_direction, count_graph));

    try_restart = !try_restart;

    if (files_sys::bin::WriteSimple(glb_files.graph_address + std::to_string(cur_direction) + ".bin", graph)) {
      RETURN_ERR("file_graph is not opened for writing\n");
    }

    WRITE_LOG("Id_proc: %d. Graph construction in the direction %d is completed, t= %lf c. \n", myid, cur_direction,
              t + omp_get_wtime());
  }

  bound_trace.clear();

  t += omp_get_wtime();
  WRITE_LOG("Full time: %lf\n", t);

  return e_completion_success;
}

#endif //! BUILD_GRAPH