/**
 * @file trace_to_face.h
 * @brief Перенос трассировки с ячеек на грани
 *
 */

#if !defined TRACE_STRUCT_H && defined MAKE_TRACE
#define TRACE_STRUCT_H

#include "solvers_struct.h"

namespace trace {

/**
 * @brief Формирует массив номеров граней, по которым определяется граничное излучение
 *
 * @param[in] grid сетка со связью граней с ячейками
 * @param[in] neighbours массив соседей
 * @param[in] face_states признаки входящих/выходящих граней
 * @param[in] sorted_id_cell упорядоченные номера ячеек
 * @param[out] graph_bound_faces номера граней
 */
void GetBoundGraphFaces(const grid_t &grid, const std::vector<IntId> &neighbours,
                        const std::vector<bits_flag_t> &face_states, const std::vector<IntId> &sorted_id_cell,
                        std::vector<IntId> &graph_bound_faces);

/**
 * @brief Get the Graph Faces object
 *
 * @param grid сетка со связью граней с ячейками
 * @param face_states признаки входящих/выходящих граней
 * @param sorted_id_cell упорядоченные номера ячеек
 * @param graph_cell_faces упорядоченные пары ячейка-граней
 */
void GetGraphFaces(const grid_t &grid,
                   const std::vector<bits_flag_t> &face_states, const std::vector<IntId> &sorted_id_cell,
                   std::vector<IntId> &graph_cell_faces);
} // namespace trace
#endif //! TRACE_STRUCT_H