#if 0 // def BUILD_GRAPH

#include <omp.h>

#include "dbgdef.h"
#include "graph_calc_omp.h"
#include "graph_config.h"
#include "graph_init_state.h"
#include "graph_struct.h"
#include "intersections.h"
#include "reader_bin.h"
#include "reader_txt.h"
#include "writer_bin.h"

namespace graph_omp {

#ifdef GRID_WITH_INNER_BOUNDARY

namespace in = intersection;
using namespace graph;
//-------------------------------------------------------------------------------------

static boundary_trace_t buf_try;

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
    in::IntersectionWithPlane(plane_face.face, vertex1, direction, buf_try.x[0]);
    if (in::InTriangle(plane_face.face, normals[in_id].n[plane_face.face_id%SIZE_CELL], buf_try.x[0]))
    // if (in_id != cur_cell && ((intersect_point - vertex1).dot(direction) < 0))
    {
      /*	std::bitset<4>id;
              FindInAndOutFaces(direction, in_id, normals, id);
              if (id[plane_face.face_id % 4] == 0) break;*/

      buf_try.id[0] = plane_face.face_id;
      id[0] = in_id;
      break;
    }
  }

  if (id[0] == -1)
    return 1;

  for (auto &in_id : inner_bound) {
    plane_face = inner_cells.find(in_id)->second;
    in::IntersectionWithPlane(plane_face.face, vertex2, direction, buf_try.x[1]);
    if (in::InTriangle(plane_face.face, normals[in_id].n[plane_face.face_id%SIZE_CELL], buf_try.x[1]))
    // if (in_id != cur_cell && ((intersect_point - vertex2).dot(direction) < 0))
    {
      /*	std::bitset<4>id;
              FindInAndOutFaces(direction, in_id, normals, id);
              if (id[plane_face.face_id % 4] == 0) break;*/
      id[1] = in_id;
      buf_try.id[1] = plane_face.face_id;
      break;
    }
  }

  if (id[1] == -1)
    return 1;

  for (auto &in_id : inner_bound) {
    plane_face = inner_cells.find(in_id)->second;
    in::IntersectionWithPlane(plane_face.face, vertex3, direction, buf_try.x[2]);
    if (in::InTriangle(plane_face.face, normals[in_id].n[plane_face.face_id%SIZE_CELL], buf_try.x[2]))
    // if (in_id != cur_cell && ((intersect_point - vertex3).dot(direction) < 0))
    {
      /*	std::bitset<4>id;
              FindInAndOutFaces(direction, in_id, normals, id);
              if (id[plane_face.face_id % 4] == 0) break;*/
      id[2] = in_id;
      buf_try.id[2] = plane_face.face_id;
      break;
    }
  }

  if (id[2] == -1) {
    return 1;
  }

  buf_try.s[0] = (vertex1 - buf_try.x[0]).norm();
  buf_try.s[1] = (vertex2 - buf_try.x[1]).norm();
  buf_try.s[2] = (vertex3 - buf_try.x[2]).norm();

  return 0;
}

static int FindCurCellWithHole(const std::vector<IntId> &next_step_el,
                               const std::vector<State> &count_in_face, std::vector<State> &count_knew_face,
                               std::list<IntId> &cur_el,
                               const std::set<IntId> &inner_part, std::set<IntId> &outter_part,
                               const std::map<IntId, FaceCell> &inner_cells,
                               const std::vector<IntId> &all_pairs_id, const Vector3 &direction,
                               const std::vector<Normals> &normals) {

  cur_el.clear();

  const int N = next_step_el.size();
#pragma omp parallel default(none) shared(next_step_el, count_in_face, count_knew_face, cur_el, inner_part, outter_part, all_pairs_id, direction, normals, inner_cells, N)
  {
    int cell;
#pragma omp for
    for (int i = 0; i < N; ++i) {
      cell = next_step_el[i];
      if (outter_part.count(cell) != 0) {
        if (count_in_face[cell] == count_knew_face[cell] + 1) { // граница не определена
          IntId try_id[3] = {-1, -1, -1};

          FindIdCellInBoundary(direction, inner_part, inner_cells, normals, cell, try_id);

          if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1)
            int a;
          else if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
                   count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
                   count_in_face[try_id[2]] == count_knew_face[try_id[2]]) { // если грань на другом конце определена
#pragma omp critical
            {
              if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
                  count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
                  count_in_face[try_id[2]] == count_knew_face[try_id[2]]) {
                cur_el.push_back(cell);
                outter_part.erase(cell);
                count_knew_face[cell]++;
              }
            }
          }
        } else if (count_in_face[cell] == count_knew_face[cell]) { // граница определена
#pragma omp critical
          {
            if (count_in_face[cell] == count_knew_face[cell]) {
              cur_el.push_back(cell);
              outter_part.erase(cell);
            }
          }
        }
      } else if (count_in_face[cell] == count_knew_face[cell]) {
#pragma omp critical
        {
          if (count_in_face[cell] == count_knew_face[cell])
            cur_el.push_back(cell);
        }
      }
    }
  }

  if (cur_el.size() == 0) {

    // плохо, но как есть. Если не смогли найти ни одну ячейку кандидата
    //попробовать пройти отдельно по внутренней границе,

    std::list<IntId> buf_erase;

    for (auto cell : outter_part) {
      if (count_in_face[cell] == count_knew_face[cell] + 1) { // граница не определена
        IntId try_id[3] = {-1, -1, -1};
        FindIdCellInBoundary(direction, inner_part, inner_cells, normals, cell, try_id);
        if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1)
          continue;
        if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
            count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
            count_in_face[try_id[2]] == count_knew_face[try_id[2]]) { // если грань на другом конце определена

          cur_el.push_back(cell);
          count_knew_face[cell]++;
          buf_erase.push_back(cell); // outter_part.erase(cell);
          continue;
        }
      } else if (count_in_face[cell] == count_knew_face[cell]) {
        buf_erase.push_back(cell); // outter_part.erase(cell);
        cur_el.push_back(cell);    //добавляет все найденные ранее границы!!!
      }
    }

    for (auto el : buf_erase)
      outter_part.erase(el);

    if (cur_el.size() == 0) {
      std::cout << "NextCell is -1\n";
      return -1;
    }
  }

  return 0;
}
#else
static int FindCurCell(const std::vector<IntId> &next_step_el,
                       const std::vector<State> &count_in_face,
                       const std::vector<State> &count_knew_face,
                       std::list<IntId> &cur_el) {

  cur_el.clear();

  const int N = next_step_el.size();
#pragma omp parallel default(none) shared(next_step_el, count_in_face, count_knew_face, cur_el, N)
  {
#pragma omp for
    for (int i = 0; i < N; i++) {
      int cell = next_step_el[i];
      if (count_in_face[cell] == count_knew_face[cell]) {
#pragma omp critical
        {
          if (count_in_face[cell] == count_knew_face[cell])
            cur_el.push_back(cell);
        }
      }
    }
  }

  if (cur_el.size() == 0) {
    std::cout << "NextCell is -1\n";
    return -1;
  }

  return 0;
}
#endif

// число входящих граней для каждой ячейки + число известных из них + начальная граница
static int FindNumberOfAllInAndDefFaces(const Vector3 &dir, const std::vector<Normals> &normals,
                                        const std::vector<State> &faces_state,
                                        std::vector<State> &count_in_face, std::vector<State> &count_knew_face,
                                        std::vector<IntId> &next_step_el) {

  const int N = normals.size(); // число ячеек

  count_in_face.assign(N, 0);
  count_knew_face.assign(N, 0);
  next_step_el.clear();

  uint32_t state;
  for (int i = 0; i < N; i++) {
    in::FindInAndOutFaces(dir, normals[i], state);
    for (int j = 0; j < 4; j++) {
      if (CHECK_BIT(state, j)) { // для i-ой ячейки добавляем:
        count_in_face[i]++;      // входную грань
        if (faces_state[i * 4 + j]) {
          count_knew_face[i]++;      //определнную входную грань (т.е. граничную)
          next_step_el.push_back(i); // начальный набор кандидатов
        }
      }
    }
  }
  return 0;
}

static int NewStep(const std::vector<IntId> &all_pairs_id,
                   const std::vector<State> &count_in_face,
                   std::vector<State> &count_knew_face,
                   std::list<IntId> &cur_el,
                   std::vector<IntId> &next_step_el) {

  int buf_size = next_step_el.size();
  next_step_el.clear();
  next_step_el.reserve(buf_size); // резервируем память на основе предыдущего шага(предпологая, что порядок величины будет тот же)

  int N = cur_el.size();

  std::set<IntId> next_step;

  for (auto cell : cur_el) {
    // по всем соседям
    for (size_t j = 0; j < 4; j++) {
      int neighbour = all_pairs_id[cell * 4 + j];
      if (neighbour < 0)
        continue;
      neighbour /= 4;

      if (count_in_face[neighbour] > count_knew_face[neighbour]) {
        count_knew_face[neighbour]++; // всегда ли эта грань будет входящей (проверка по нормалям??)
        if (next_step.count(neighbour) == 0) {
          next_step.emplace(neighbour);
          next_step_el.push_back(neighbour); // ячейка была изменена, проверить ее готовность на следующем шаге
        }
      }
    }
  }

  return 0;
}

int RunGraphModule() {
  std::vector<IntId> neighbours;          ///< соседние ячейки
  std::set<IntId> inter_boundary_face_id; ///< id внутренних граней [i * CELL_SIZE + j]
  std::vector<Normals> normals;           ///< нормали
  std::map<IntId, FaceCell> inter_faces;  ///< внутренние грани с ключом-номером ячейки
  grid_directions_t grid_dir;

  uint32_t err = 0;
  err |= files_sys::bin::ReadSimple(glb_files.name_file_neigh, neighbours);
  err |= files_sys::txt::ReadSimple(glb_files.base_address + F_INTERNAL_BOUND, inter_boundary_face_id);
  err |= files_sys::bin::ReadNormals(glb_files.base_address + F_NORMALS, normals);
  err |= files_sys::txt::ReadInitBoundarySetInFaces(glb_files.base_address + F_FACE_ID, inter_faces);
  err |= files_sys::txt::ReadSphereDirectionСartesian(glb_files.name_file_sphere_direction, grid_dir);

  if (err != 0) {
    RETURN_ERR("error during reading\n");
  }

  std::vector<State> count_in_face(normals.size(), 0);  ///< число входящих граней ячейки
  std::vector<State> count_def_face(normals.size(), 0); ///< число определённых граней ячейки
  std::vector<IntId> graph(normals.size(), 0);          ///< упорядоченный набор ячеек

  std::set<IntId> inner_part; ///< часть граничных ячеек через которые луч покидай расчётную область
  std::set<IntId> outer_part; ///< часть граничных ячеек через которые луч возвращается в расчётную область

  std::vector<State> faces_state;

  int myid = get_mpi_id();
  int np = get_mpi_np();
  bool flag = true;
  std::list<IntId> cur_el_OMP;         ///<текущая границы
  std::vector<IntId> next_step_el_OMP; ///< кандидаты на следующую границу

  double t = -omp_get_wtime();
  for (int cur_direction = myid; cur_direction < grid_dir.size; cur_direction += np) {

    WRITE_LOG("Direction #: %d\n", cur_direction);

    InitFacesState(neighbours, inter_faces, faces_state);
    Vector3 direction = grid_dir.directions[cur_direction].dir;

    DivideInnerBoundary(direction, normals, inter_boundary_face_id, inner_part, outer_part);

    bound_trace.clear();
    bound_trace.reserve(outer_part.size() * 3); // потенциально в массив могут войдут все ячейки внутренней границы

    FindNumberOfAllInAndDefFaces(direction, normals, faces_state, count_in_face, count_def_face,
                                 next_step_el_OMP);

    int count_graph = 0;              // число ячеек вошедших в граф
    graph.assign(normals.size(), -1); // для отлавливания ошибочных направлений
    bool try_restart = true;

    while (next_step_el_OMP.size() && flag) {

#ifdef GRID_WITH_INNER_BOUNDARY
      IntId cur_ret = FindCurCellWithHole(next_step_el_OMP, count_in_face, count_def_face, cur_el_OMP,
                                          inner_part, outer_part, inter_faces, neighbours, direction, normals);
#else
      IntId cur_ret = FindCurCell(next_step_el_OMP, count_in_face, count_def_face, cur_el_OMP);
#endif // GRID_WITH_INNER_BOUNDARY

      if (cur_ret == -1) {
        WRITE_LOG("Warning proc: %d, dir= %d, processed %d cells", myid, cur_direction, count_graph);

        cur_el_OMP.clear();
        std::vector<IntId> foo;
        std::set<IntId> next;

        if (TryRestart(count_in_face, count_def_face, outer_part, foo, next) == e_completion_success) {
          WRITE_LOG("\n\n Warning!!! try_restart %d \n\n", cur_direction);

          next_step_el_OMP.assign(next.begin(), next.end());
          continue;
        }
      }

      NewStep(neighbours, count_in_face, count_def_face, cur_el_OMP, next_step_el_OMP);

      for (auto el : cur_el_OMP) {
        graph[count_graph] = el;
        count_graph++;
      }

    } // while

    DIE_IF_ACTION(count_graph < graph.size(), WRITE_LOG_ERR("Error size graph[%d] %d\n", cur_direction, count_graph));

    try_restart = !try_restart;

    if (files_sys::bin::WriteSimple(glb_files.graph_address + F_GRAPH + std::to_string(cur_direction) + ".bin", graph)) {
      RETURN_ERR("file_graph is not opened for writing\n");
    }

    WRITE_LOG("Id_proc: %d. Graph construction in the direction %d is completed, t= %lf c. \n", myid, cur_direction,
              t + omp_get_wtime());
  }
}

} // namespace graph_omp

#endif //! BUILD_GRAPH