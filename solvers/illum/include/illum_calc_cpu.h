/**
 * @file illum_calc_cpu.h
 * @brief Функции расчета излучения без видеокарты с использованием одного узла и нескольких потоков
 *
 */
#if !defined ILLUM_PART_H && defined ILLUM && defined SOLVERS
#define ILLUM_PART_H
#include "geo_types.h"
#include "solvers_struct.h"

namespace illum {

/**
 * @brief Пространство имён расчёта излучения на cpu  без видеокарты
 *
 */
namespace cpu {

int CalculateIllum(const grid_directions_t &grid_direction, const std::vector<std::vector<bits_flag_t>> &face_states,
                   const std::vector<IntId> &neighbours, const std::vector<std::vector<IntId>> &inner_bound_code,
                   const std::vector<std::vector<cell_local>> &vec_x0, std::vector<BasePointTetra> &vec_x,
                   const std::vector<std::vector<IntId>> &sorted_id_cell, grid_t &grid);

void CalculateIllumParam(const grid_directions_t &grid_direction, grid_t &grid);

int CalculateAdditionalIllum(const grid_directions_t &grid_direction, const std::vector<std::vector<bits_flag_t>> &face_states,
                             const std::vector<IntId> &neighbours, const std::vector<std::vector<IntId>> &inner_bound_code,
                             const std::vector<std::vector<cell_local>> &vec_x0, std::vector<BasePointTetra> &vec_x,
                             const std::vector<std::vector<IntId>> &sorted_id_cell, grid_t &grid);

int CalculateIllumFace(const grid_directions_t &grid_direction,
                       const std::vector<std::vector<IntId>> &inner_bound_code,
                       const std::vector<std::vector<cell_local>> &vec_x0,
                       const std::vector<std::vector<graph_pair_t>> &sorted_graph,
                       const std::vector<std::vector<IntId>> &sorted_id_bound_face,
                       grid_t &grid);

} // namespace cpu

} // namespace illum

#endif //! ILLUM_PART_H