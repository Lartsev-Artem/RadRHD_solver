#if !defined GRAPH_STRUCT_H && defined BUILD_GRAPH
#define GRAPH_STRUCT_H

#include <Eigen/Dense>
#include <vector>

/*! \addtogroup graph Модуль построения графов
    \brief Модуль содержит функции построения графов для маршевого алгоритма расчёта излучения
    @{
*/

/**
 * @brief Пространство имён модуля построения графов
 *
 */
namespace graph {

extern std::vector<int> id_try_surface;            // id граней, определяющих внутренюю границу
extern std::vector<double> dist_try_surface;       // расстояния между точками (через полость внутри)
extern std::vector<Eigen::Vector3d> x_try_surface; // x точка выхода

extern uint64_t id_try_size;
extern uint64_t dist_try_size;
extern uint64_t x_try_size;

struct try_solve_t {
  int id_1;
  int id_2;
  int id_3;

  double s_1;
  double s_2;
  double s_3;

  Eigen::Vector3d x1;
  Eigen::Vector3d x2;
  Eigen::Vector3d x3;
};
extern try_solve_t buf_try;

/**
 * @brief код состояние грани
 *
 */
enum e_face_state_t {
  e_face_state_undef = 0, ///< грань требует определения через соседнии
  e_face_state_def = 1    ///< грань определена через соседнии
};

} // namespace graph
#endif //! GRAPH_STRUCT_H
