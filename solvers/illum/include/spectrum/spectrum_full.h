#if !defined SPECTRUM_FULL_UTILS_H && defined SPECTRUM
#define SPECTRUM_FULL_UTILS_H

#include "solvers_struct.h"
#ifdef SAVE_FULL_SPECTRUM

namespace illum {
namespace full_spectrum {

Type BoundaryConditions(const IdType type_bound, const Type frq0, const Type frq);

Type ReCalcIllum(const IdType num_dir, const std::vector<std::vector<Type>> &inter_coef, grid_t &grid, const IdType dir_disp);

int CalculateIllum(const grid_directions_t &grid_direction,
                   const std::vector<std::vector<IntId>> &inner_bound_code,
                   const std::vector<align_cell_local> &vec_x0,
                   const std::vector<std::vector<graph_pair_t>> &sorted_graph,
                   const std::vector<std::vector<IntId>> &sorted_id_bound_face,
                   grid_t &grid);
} // namespace full_spectrum
} // namespace illum

#endif //! SAVE_FULL_SPECTRUM
#endif //! SPECTRUM_FULL_UTILS_H &&  SPECTRUM