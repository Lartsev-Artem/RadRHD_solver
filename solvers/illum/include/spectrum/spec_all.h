#ifndef SPEC_ALL_H
#define SPEC_ALL_H
#include "geo_types.h"
#include "solvers_struct.h"

/*! \addtogroup illum Модуль расчёта излучения
    @{
*/

namespace illum {

/**
 * @brief Пространство имён расчёта излучения на cpu  без видеокарты
 *
 */
namespace spec {

int CalculateIllumFace(const grid_directions_t &grid_direction,
                       const std::vector<std::vector<IntId>> &inner_bound_code,
                       const std::vector<align_cell_local> &vec_x0,
                       const std::vector<std::vector<graph_pair_t>> &sorted_graph,
                       const std::vector<std::vector<IntId>> &sorted_id_bound_face,
                       grid_t &grid);

Type BoundaryConditions(const IdType type_bound, const Type frq0, const Type frq1);
Type GetIllum(const Vector3 &dir, const Vector3 &x,
              const Type s,
              const Type I_0,
              const Type int_scattering,
              const Type frq, const Type frq0,
              elem_t &cell);

Type ReCalcIllum(const IdType num_dir, const IdType num_frq, const std::vector<Type> &inter_coef, grid_t &grid, IdType mpi_dir_shift = 0);

int RunIllumSpectrumModule();
} // namespace spec

} // namespace illum
#endif //! SPEC_ALL_H