#ifdef BUILD_GRAPH
#include "graph_calc.h"

#include "graph_config.h"
#include "graph_struct.h"
#include "intersections.h"

#ifdef GRID_WITH_INNER_BOUNDARY
namespace in = intersection;

/**
 * @brief Функция ищет пересечения луча из выходящей части внутренней границы с входящей
 *
 * @param[in] out_cell номер ячейки на выходящей части
 * @param[in] direction направление
 * @param[in] inner_part часть определяющей границы
 * @param[in] inner_faces все граничные грани(внутренней границы)
 * @param[in] normals нормали
 * @param[out] trace параметры луча
 * @return int число найденных пересечений
 * @warning не проверят известна ли inner_part на данный момент. Проверяет только геометрию
 */
static int FindCellOnInnerPartBoundary(const int out_cell,
                                       const Vector3 &direction,
                                       const std::set<IntId> &inner_part,
                                       const std::map<IntId, FaceCell> &inner_faces,
                                       const std::vector<Normals> &normals,
                                       graph::boundary_trace_t &trace) {

  // вершины треугольника.
  Face face = inner_faces.find(out_cell)->second.face;
  Vector3 P1(face.A.data());
  Vector3 P2(face.B.data());
  Vector3 P3(face.C.data());

  // середины сторон (противолежащий точке P1,P2,P3 соответственно)
  Vector3 P11 = (P3 + P2) / 2;
  Vector3 P22 = (P3 + P1) / 2;
  Vector3 P33 = (P2 + P1) / 2;

  // точки на медианах
  Vector3 vertex[3] = {P1 + (P11 - P1) / 3,
                       P2 + (P22 - P2) / 3,
                       P3 + (P33 - P3) / 3};

  int count_intersection = 0;

  /// \note обратная трассировка. Ищем пересечение от точки вхождения до выхождения
  for (int i = 0; i < 3; i++) {
    for (auto &in_cell : inner_part) {

      FaceCell start_face = inner_faces.find(in_cell)->second; //грань на противолежащей стороне

      in::IntersectionWithPlane(start_face.face, vertex[i], -direction, trace.x[i]);
      if (in::InTriangle(start_face.face_id, start_face.face, normals[in_cell], trace.x[i])) {
        trace.id[i] = start_face.face_id;
        trace.s[i] = (vertex[i] - trace.x[i]).norm();
        count_intersection++;
        break;
      }
    }
  }
  return count_intersection;
}

int graph::FindCurFrontWithHole(const Vector3 &direction,
                                const std::vector<Normals> &normals,
                                const std::map<IntId, FaceCell> &inner_faces,
                                const std::set<IntId> &next_candidates,
                                const std::vector<State> &count_in_face,
                                const std::set<IntId> &inner_part,
                                std::vector<IntId> &cur_front,
                                std::vector<State> &count_def_face,
                                std::set<IntId> &outer_part) {

  cur_front.clear(); // очищаем текущую границу

  for (auto cell_id : next_candidates) {

    // ячейка на границе выхода луча из подобласти
    if (outer_part.count(cell_id) != 0) {

      // граница не определена, но осталась последняя граничная часть
      if (count_in_face[cell_id] == count_def_face[cell_id] + 1) {

        boundary_trace_t trace;

        if (FindCellOnInnerPartBoundary(cell_id, direction, inner_part, inner_faces, normals, trace) != 3) {
          D_LD; //!не нашли определяющей геометрии
        }

        int count = 0;
        for (int i = 0; i < 3; i++) {
          int cell_id = trace.id[i] / CELL_SIZE;
          count += (count_in_face[cell_id] == count_def_face[cell_id]);
        }

        // если грань на другом конце полностью определена
        if (count == 3) {

          bound_trace.push_back(trace);
          cur_front.push_back(cell_id);
          outer_part.erase(cell_id);
          count_def_face[cell_id]++;
          continue;
        }
      }
      //! не должно быть так.  // граница определена
      // else if (count_in_face[cell_id] == count_def_face[cell_id]) {
      //   cur_front.push_back(cell_id);
      //   outer_part.erase(cell_id);
      //}
    } else if (count_in_face[cell_id] == count_def_face[cell_id]) {
      cur_front.push_back(cell_id);
    }
  }

  if (cur_front.size() == 0)
    return e_completion_fail;
  return e_completion_success;
}

#else
int graph::FindCurFront(const std::set<IntId> &next_candidates,
                        const std::vector<State> &count_in_face,
                        const std::vector<State> &count_def_face,
                        std::vector<IntId> &cur_front) {

  cur_front.clear(); // очищаем текущую границу

  for (auto cell : next_candidates) {
    if (count_in_face[cell] == count_def_face[cell])
      cur_front.push_back(cell);
  }

  if (cur_front.size() == 0) {
    return e_completion_fail;
  }

  return e_completion_success;
}

#endif
void graph::NewStep(const std::vector<IntId> &neighbours,
                    const std::vector<State> &count_in_face,
                    const std::vector<IntId> &cur_front,
                    std::vector<State> &count_def_face,
                    std::set<IntId> &next_candidate) {

  next_candidate.clear(); //очищаем список кандидатов

  for (auto cell : cur_front) {
    // по всем соседям
    for (int j = 0; j < CELL_SIZE; j++) {
      int neigh = neighbours[cell * CELL_SIZE + j];
      if (neigh < 0)
        continue;
      neigh /= CELL_SIZE;

      // ячейка была изменена, проверить ее готовность на следующем шаге
      if (count_in_face[neigh] > count_def_face[neigh]) {
        next_candidate.emplace(neigh);
        count_def_face[neigh]++;
      }
    }
  }
}

#endif //! BUILD_GRAPH