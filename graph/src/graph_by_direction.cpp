#if defined BUILD_GRAPH
#include "graph_config.h"
#include "graph_main.h"


#include "dbgdef.h"

#include "graph_calc.h"
#include "graph_init_state.h"
#include "graph_inner_bound.h"
#include "graph_struct.h"

int graph::RunGraphModule(TracerData &data, const Vector3 &direction) {
  const std::vector<IntId> &neighbours = data.neighbours;
  const std::set<IntId> &inter_boundary_face_id = data.inter_boundary_face_id;
  const std::vector<Normals> &normals = data.normals;
  const std::map<IntId, FaceCell> &inter_faces = data.inter_faces;
  std::vector<IntId> &graph = data.graph;

  const size_t num_cells = normals.size();

  std::vector<State> count_in_face(num_cells, 0);  ///< число входящих граней ячейки
  std::vector<State> count_def_face(num_cells, 0); ///< число определённых граней ячейки

  std::set<IntId> inner_part; ///< часть граничных ячеек через которые луч покидай расчётную область
  std::set<IntId> outer_part; ///< часть граничных ячеек через которые луч возвращается в расчётную область

  std::vector<State> faces_state; ///< состояние граней (определена не определена)

  std::vector<IntId> cur_el;    ///<текущая границы
  std::set<IntId> next_step_el; ///< кандидаты на следующую границу

  for (int cur_direction = 0; cur_direction < 1; cur_direction++) {
    InitFacesState(neighbours, inter_faces, faces_state);

    DivideInnerBoundary(direction, normals, inter_boundary_face_id, inner_part, outer_part);

    bound_trace.clear();
    bound_trace.reserve(outer_part.size() * 3); // потенциально в массив могут войдут все ячейки внутренней границы

    FindNumberOfAllInAndDefFaces(direction, normals, faces_state, count_in_face, count_def_face, next_step_el);

    int count_graph = 0;         // число ячеек вошедших в граф
    graph.assign(num_cells, -1); // для отлавливания ошибочных направлений
    bool try_restart = true;

    while (next_step_el.size()) {

#ifdef GRID_WITH_INNER_BOUNDARY
      int cur_ret = FindCurFrontWithHole(direction, normals, inter_faces, next_step_el,
                                         count_in_face, inner_part, cur_el, count_def_face, outer_part);
#else
      int cur_ret = FindCurFront(next_step_el, count_in_face, count_def_face, cur_el);
#endif // GRID_WITH_INNER_BOUNDARY

      NewStep(neighbours, count_in_face, cur_el, count_def_face, next_step_el);

      for (auto el : cur_el) {
        graph[count_graph] = el;
        count_graph++;
      }

      if ((graph.size() != count_graph && next_step_el.size() == 0) || cur_ret != e_completion_success) {
        WRITE_LOG("Warning proc: dir= %d, processed %d cells\n", cur_direction, count_graph);

        if (TryRestart(count_in_face, count_def_face, outer_part, cur_el, next_step_el) == e_completion_success) {
          // WRITE_LOG("\n\n Warning!!! try_restart %d \n\n", cur_direction);
          continue;
        }
      }

    } // while

    DIE_IF_ACTION(count_graph < graph.size(), WRITE_LOG_ERR("Error size graph[%d] %d\n", cur_direction, count_graph));
  }

  bound_trace.clear();

#if defined GRID_WITH_INNER_BOUNDARY && defined USE_CUDA && defined GRAPH_TRACING_INNER_BOUNDARY
  trace_through_boundary::ClearDevice();
#endif

  return e_completion_success;
}

#endif //! BUILD_GRAPH