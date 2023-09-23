#ifdef BUILD_GRAPH
#include "graph_calculation.h"
#include "graph_struct.h"
#include "intersections.h"

#include <list>

namespace graph {

#ifdef GRID_WITH_INNER_BOUNDARY
namespace in = intersection;
static int FindIdCellInBoundary(const Vector3 &direction,
                                const std::set<IntId> &inner_bound,
                                const std::map<IntId, FaceCell> &inner_cells,
                                const std::vector<Normals> &normals,
                                const int cur_cell, int *id) {

  // вершины треугольника.
  Face face = inner_cells.find(cur_cell)->second.face;
  Vector3 P1(face.A.data());
  Vector3 P2(face.B.data());
  Vector3 P3(face.C.data());

  // середины сторон (противолежащий точке P1,P2,P3 соответственно)
  Vector3 P11 = (P3 + P2) / 2;
  Vector3 P22 = (P3 + P1) / 2;
  Vector3 P33 = (P2 + P1) / 2;

  // точки на медианах
  Vector3 vertex1 = P1 + (P11 - P1) / 3;
  Vector3 vertex2 = P2 + (P22 - P2) / 3;
  Vector3 vertex3 = P3 + (P33 - P3) / 3;

  // ищем пересечения vertex->direction c гранями внутренней границы
  // Vector3 intersect_point1;  //buf_try.x1
  // Vector3 intersect_point2;  //buf_try.x2
  // Vector3 intersect_point3;  //buf_try.x3

  FaceCell plane_face;
  for (auto &in_id : inner_bound) {
    plane_face = inner_cells.find(in_id)->second;
    in::IntersectionWithPlane(plane_face.face, vertex1, direction, buf_try.x1);
    if (in::InTriangle(plane_face.face_id, plane_face.face, normals[in_id], buf_try.x1))
    // if (in_id != cur_cell && ((intersect_point - vertex1).dot(direction) < 0))
    {
      /*	std::bitset<4>id;
              FindInAndOutFaces(direction, in_id, normals, id);
              if (id[plane_face.face_id % 4] == 0) break;*/

      buf_try.id_1 = plane_face.face_id;
      id[0] = in_id;
      break;
    }
  }

  if (id[0] == -1)
    return 1;

  for (auto &in_id : inner_bound) {
    plane_face = inner_cells.find(in_id)->second;
    in::IntersectionWithPlane(plane_face.face, vertex2, direction, buf_try.x2);
    if (in::InTriangle(plane_face.face_id, plane_face.face, normals[in_id], buf_try.x2))
    // if (in_id != cur_cell && ((intersect_point - vertex2).dot(direction) < 0))
    {
      /*	std::bitset<4>id;
              FindInAndOutFaces(direction, in_id, normals, id);
              if (id[plane_face.face_id % 4] == 0) break;*/
      id[1] = in_id;
      buf_try.id_2 = plane_face.face_id;
      break;
    }
  }

  if (id[1] == -1)
    return 1;

  for (auto &in_id : inner_bound) {
    plane_face = inner_cells.find(in_id)->second;
    in::IntersectionWithPlane(plane_face.face, vertex3, direction, buf_try.x3);
    if (in::InTriangle(plane_face.face_id, plane_face.face, normals[in_id], buf_try.x3))
    // if (in_id != cur_cell && ((intersect_point - vertex3).dot(direction) < 0))
    {
      /*	std::bitset<4>id;
              FindInAndOutFaces(direction, in_id, normals, id);
              if (id[plane_face.face_id % 4] == 0) break;*/
      id[2] = in_id;
      buf_try.id_3 = plane_face.face_id;
      break;
    }
  }

  if (id[2] == -1)
    return 1;

  buf_try.s_1 = (vertex1 - buf_try.x1).norm();
  buf_try.s_2 = (vertex2 - buf_try.x2).norm();
  buf_try.s_3 = (vertex3 - buf_try.x3).norm();

  return 0;
}

int FindCurCellWithHole(const Vector3 &direction,
                        const std::vector<Normals> &normals,
                        const std::map<IntId, FaceCell> &inner_faces,
                        const std::set<IntId> &next_step_el,
                        const std::vector<State> &count_in_face,
                        const std::set<IntId> &inner_part,
                        std::vector<IntId> &cur_front,
                        std::vector<State> &count_def_face,
                        std::set<IntId> &outer_part) {

  cur_front.clear(); // очищаем текущую границу

  for (auto cell : next_step_el) {
    if (outer_part.count(cell) != 0) {

      // граница не определена
      if (count_in_face[cell] == count_def_face[cell] + 1) {

        IntId try_id[3] = {-1, -1, -1};

        if (FindIdCellInBoundary(direction, inner_part, inner_faces, normals, cell, try_id))
          continue;

        if (count_in_face[try_id[0]] == count_def_face[try_id[0]] &&
            count_in_face[try_id[1]] == count_def_face[try_id[1]] &&
            count_in_face[try_id[2]] == count_def_face[try_id[2]]) { // если грань на другом конце определена

          // id_try_surface.push_back(buf_try.id_1);
          // id_try_surface.push_back(buf_try.id_2);
          // id_try_surface.push_back(buf_try.id_3);

          // x_try_surface.push_back(buf_try.x1);
          // x_try_surface.push_back(buf_try.x2);
          // x_try_surface.push_back(buf_try.x3);

          // dist_try_surface.push_back(buf_try.s_1);
          // dist_try_surface.push_back(buf_try.s_2);
          // dist_try_surface.push_back(buf_try.s_3);

          cur_front.push_back(cell);
          outer_part.erase(cell);
          count_def_face[cell]++;
          continue;
        }
      } else if (count_in_face[cell] == count_def_face[cell]) { // граница определена
        cur_front.push_back(cell);
        outer_part.erase(cell);
      }
    } else if (count_in_face[cell] == count_def_face[cell]) {
      cur_front.push_back(cell);
    }
  }

  if (cur_front.size() == 0) {

    // плохо, но как есть. Если не смогли найти ни одну ячейку кандидата \
		попробовать пройти отдельно по внутренней границе,

    std::list<IntId> buf_erase;

    for (auto cell : outer_part) {
      if (count_in_face[cell] == count_def_face[cell] + 1) { // граница не определена
        IntId try_id[3] = {-1, -1, -1};
        FindIdCellInBoundary(direction, inner_part, inner_faces, normals, cell, try_id);
        if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1)
          continue;
        if (count_in_face[try_id[0]] == count_def_face[try_id[0]] &&
            count_in_face[try_id[1]] == count_def_face[try_id[1]] &&
            count_in_face[try_id[2]] == count_def_face[try_id[2]]) { // если грань на другом конце определена

          cur_front.push_back(cell);
          count_def_face[cell]++;
          buf_erase.push_back(cell); // outer_part.erase(cell);
          continue;
        }
      } else if (count_in_face[cell] == count_def_face[cell]) {
        buf_erase.push_back(cell); // outer_part.erase(cell);
        cur_front.push_back(cell);
      }
    }

    for (auto el : buf_erase)
      outer_part.erase(el);

    if (cur_front.size() == 0) {
      std::cout << "NextCell is -1\n";
      return -1;
    }
  }

  return 0;
}

#endif
int FindCurFront(const std::set<IntId> &next_step_el,
                 const std::vector<State> &count_in_face,
                 const std::vector<State> &count_def_face,
                 std::vector<IntId> &cur_front) {

  cur_front.clear();

  for (auto cell : next_step_el) {
    if (count_in_face[cell] == count_def_face[cell])
      cur_front.push_back(cell);
  }

  if (cur_front.size() == 0) {
    return e_completion_fail;
  }

  return e_completion_success;
}
void NewStep(const std::vector<IntId> &neighbours,
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

} // namespace graph
#endif //! BUILD_GRAPH