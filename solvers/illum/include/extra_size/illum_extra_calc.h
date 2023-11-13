#if !defined ILLUM_EXTRA_MAIN_H && defined ILLUM && defined SOLVERS
#define ILLUM_EXTRA_MAIN_H

#include "solvers_struct.h"
namespace illum {
namespace extra_size {

int CalculateIllum(const grid_directions_t &grid_direction, const std::vector<std::vector<State>> &face_states, const std::vector<int> &pairs,
                   const std::vector<std::vector<IntId>> &inner_bound_code,
                   const std::vector<std::vector<cell_local>> &vec_x0, std::vector<BasePointTetra> &vec_x,
                   const std::vector<std::vector<int>> &sorted_id_cell, grid_t &grid);

}
} // namespace illum

#endif //! ILLUM_EXTRA_MAIN_H