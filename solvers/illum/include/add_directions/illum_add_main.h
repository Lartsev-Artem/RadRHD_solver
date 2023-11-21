#if !defined ILLUM_ADD_MAIN_H
#define ILLUM_ADD_MAIN_H

#include "geo_types.h"
#include "solvers_config.h"

namespace illum {

namespace additional_direction {
#ifndef TRANSFER_CELL_TO_FACE
int RunModule();
#endif
} // namespace additional_direction
} // namespace illum

#endif //! ILLUM_ADD_MAIN_H