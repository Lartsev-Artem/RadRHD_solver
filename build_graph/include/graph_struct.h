/**
 * @file graph_struct.h
 * @brief Объявления внутренних структур модуля
 */
#if !defined GRAPH_STRUCT_H && defined BUILD_GRAPH
#define GRAPH_STRUCT_H

#include <Eigen/Dense>
#include <vector>

/*! \addtogroup graph Модуль построения графов
    @{
*/

namespace graph {

struct boundary_trace_t {
  int id[3];            ///< номер грани в глобальной индексации, определяющей i-ый узел граничной грани
  double s[3];          ///< расстояние пройденное лучём от определяющей точки
  Eigen::Vector3d x[3]; ///< точка определяющей грани (в глобальных координатах)
};
extern std::vector<boundary_trace_t> bound_trace;

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
