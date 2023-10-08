#if !defined RHLLC_INIT_H && defined SOLVERS
#define RHLLC_INIT_H

#include "solvers_struct.h"

/*! \addtogroup rhllc Модуль расчета газовой динамики в релятивистской постановке
    @{
*/

namespace rhllc {
void SetCfgDefault(hllc_value_t &hllc_set);
int Init(std::string &file_init_value, std::vector<elem_t> &cells);
} // namespace rhllc
#endif //! RHLLC_INIT_H