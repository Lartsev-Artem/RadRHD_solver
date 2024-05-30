#if !defined SPECTRUM_CALC_H && defined SPECTRUM
#define SPECTRUM_CALC_H
#include "geo_types.h"
#include "solvers_struct.h"

/*! \addtogroup illum Модуль расчёта излучения
    @{
*/

namespace illum {

/**
 * @brief Пространство имён расчёта спектра
 *
 */
namespace spectrum {

int CalculateSpectrum(const grid_directions_t &grid_direction,
                      const std::vector<std::vector<IntId>> &inner_bound_code,
                      const std::vector<align_cell_local> &vec_x0,
                      const std::vector<std::vector<graph_pair_t>> &sorted_graph,
                      const boundary_faces_by_directions_t &boundary_faces,
                      grid_t &grid);

};

}; // namespace illum

#endif //! SPECTRUM_CALC_H