#if !defined RHLLC_UTILS_H && defined SOLVERS
#define RHLLC_UTILS_H

#include "solvers_struct.h"

/*! \addtogroup rhllc Модуль расчета газовой динамики в релятивистской постановке
    @{
*/

namespace rhllc {
void BoundConditions(const face_t &f, const std::vector<elem_t> &cells, flux_all_t &bound);
Type GetTimeStep(const hllc_value_t &hllc_set, const std::vector<elem_t> &cells);
void HllcPhysToConv(std::vector<elem_t> &cells);
} // namespace rhllc
#endif //! RHLLC_UTILS_H