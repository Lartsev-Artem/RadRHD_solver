#if defined BUILD_GRAPH
#include "graph_init_state.h"
#include "graph_struct.h"
#include "intersections.h"

#include "graph_config.h"

static uint32_t restart_counter = 0; ///< счетчик рестартов (для попытки перестроить граф не с нуля)

int graph::TryRestart(const std::vector<State> &count_in_face,
                      const std::vector<State> &count_def_face,
                      const std::set<IntId> &outer_part,
                      std::vector<IntId> &cur_front,
                      std::set<IntId> &next_candidate) {

  restart_counter++;
  cur_front.clear();
  next_candidate.clear();

#ifdef GRID_WITH_INNER_BOUNDARY
  if (outer_part.size() != 0) {

    next_candidate.insert(outer_part.begin(), outer_part.end());

    //пытаться определять внутреннюю границу (в надежде на трассировку)
    if (restart_counter < GRAPH_MAX_RESTART / 2)
      return e_completion_success;
  }
#endif

  for (int i = 0; i < count_in_face.size(); i++) {
    if (count_in_face[i] > count_def_face[i]) {
      next_candidate.emplace(i); // все неопределённые ячейки
    }
  }

  if (restart_counter > GRAPH_MAX_RESTART)
    return e_completion_fail;
  return e_completion_success;
}

void graph::InitFacesState(const std::vector<IntId> &neighbours,
                           const std::map<IntId, FaceCell> &inter_faces,
                           std::vector<State> &faces_state) {

  restart_counter = 0; /// сброс на каждом направлении

  const int n = neighbours.size();
  faces_state.assign(n, e_face_state_undef); //инициализируем состояние

  for (int i = 0; i < n; i++) {
    if (neighbours[i] < 0)
      faces_state[i] = e_face_state_def; // граница будет определена по ГУ
  }

  // сбрасываем внутреннюю границу
  for (auto &el : inter_faces)
    faces_state[el.second.face_id] = e_face_state_undef;
}

void graph::DivideInnerBoundary(const Vector3 &direction,
                                const std::vector<Normals> &normals,
                                const std::set<IntId> &inter_boundary_face_id,
                                std::set<IntId> &inner_part,
                                std::set<IntId> &outer_part) {
  inner_part.clear();
  outer_part.clear();

  for (auto id : inter_boundary_face_id) {
    bits_flag_t state = 0;

    //здесь мы определяем все грани в ячейке, но фактически нас интересует только граничная
    intersection::FindInAndOutFaces(direction, normals[id / 4], state);

    /// \note грань входящая -> граница выходящая
    /// (т.е. определяется через пересечения аккретора и диска или трассировки)

    if (CHECK_BIT(state, id % 4) == e_face_type_in)
      outer_part.emplace(id / 4);
    else
      inner_part.emplace(id / 4);
  }
}

void graph::FindNumberOfAllInAndDefFaces(const Vector3 &dir,
                                         const std::vector<Normals> &normals,
                                         const std::vector<State> &faces_state,
                                         std::vector<State> &count_in_face,
                                         std::vector<State> &count_def_face,
                                         std::set<IntId> &next_step_el) {

  const int N = normals.size(); // число ячеек

  count_in_face.assign(N, 0);
  count_def_face.assign(N, 0);
  next_step_el.clear();

  for (int i = 0; i < N; i++) {

    bits_flag_t state = 0;
    intersection::FindInAndOutFaces(dir, normals[i], state);

    for (int j = 0; j < CELL_SIZE; j++) {

      // для i-ой ячейки добавляем:
      if (CHECK_BIT(state, j) == e_face_type_in) {
        count_in_face[i]++; // входную грань
        if (faces_state[i * CELL_SIZE + j] == e_face_state_def) {
          count_def_face[i]++;     //определённую входную грань (т.е. граничную)
          next_step_el.emplace(i); // начальный набор кандидатов
        }
      }
    }
  }
}
#endif //!#if defined BUILD_GRAPH