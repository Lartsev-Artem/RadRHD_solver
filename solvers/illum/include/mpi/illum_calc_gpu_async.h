/**
 * @file illum_calc_gpu_async.h
 * @brief Вызов основного расчета излучения
 *
 */
#if !defined ILLUM_PART_H && defined ILLUM && defined SOLVERS
#define ILLUM_PART_H
#include "geo_types.h"
#include "solvers_struct.h"

/*! \addtogroup illum Модуль расчёта излучения
    @{
*/

namespace illum {

/**
 * @brief Пространство имён расчёта излучения на  видеокарте с асинхронной пересылкой между процессами
 *
 */
namespace gpu_async {

#ifndef TRANSFER_CELL_TO_FACE
int CalculateIllum(const grid_directions_t &grid_direction, const std::vector<std::vector<State>> &face_states, const std::vector<int> &pairs,
                   const std::vector<std::vector<IntId>> &inner_bound_code,
                   const std::vector<std::vector<cell_local>> &vec_x0, std::vector<BasePointTetra> &vec_x,
                   const std::vector<std::vector<int>> &sorted_id_cell, grid_t &grid);
#endif
} // namespace gpu_async

/**
 * @brief Пространство имён расчёта излучения на  видеокарте с асинхронной пересылкой между процессами и несколькими видеокартами на узле
 *
 */
namespace separate_gpu {

#ifdef TRANSFER_CELL_TO_FACE
int CalculateIllum(const grid_directions_t &grid_direction,
                   const std::vector<std::vector<IntId>> &inner_bound_code,
                   const std::vector<align_cell_local> &vec_x0,
                   const std::vector<std::vector<graph_pair_t>> &sorted_graph,
                   const boundary_faces_by_directions_t &boundary_faces,
                   grid_t &grid);

int CalculateIllumByDirection(const Vector3 &direction,
                const align_cell_local &vec_x0,
                const std::vector<graph_pair_t> &sorted_graph,
                const std::vector<IntId>& sorted_id_bound_face,
                grid_t &grid);
#endif
} // namespace separate_gpu

} // namespace illum
#endif
