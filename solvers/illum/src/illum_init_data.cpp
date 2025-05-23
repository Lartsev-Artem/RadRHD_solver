#if defined SOLVERS && defined ILLUM

#include "illum_init_data.h"

#include "global_value.h"
#include "reader_bin.h"
#include "reader_txt.h"
#include "solvers_struct.h"
int illum::InitRadiationState(const std::string &address_data, grid_t &grid) {

  std::vector<Type> data;
  if (files_sys::bin::ReadSimple(address_data + F_ABSORPCOEF, data))
    return e_completion_fail;

  if (data.size() != grid.cells.size())
    RETURN_ERR("bad size %s: %lu %lu", F_ABSORPCOEF, data.size(), grid.cells.size());

  for (size_t i = 0; i < data.size(); i++) {
    grid.cells[i].illum_val.absorp_coef = data[i];
  }

  if (files_sys::bin::ReadSimple(address_data + F_RADLOOSERATE, data))
    return e_completion_fail;

  if (data.size() != grid.cells.size())
    RETURN_ERR("bad size %s: %lu %lu", F_RADLOOSERATE, data.size(), grid.cells.size());

  for (size_t i = 0; i < data.size(); i++) {
    grid.cells[i].illum_val.rad_en_loose_rate = data[i];
  }

#if GEOMETRY_TYPE == MAIN_ELLIPSE
  if (files_sys::bin::ReadSimple(address_data + F_DENSITY, data))
    return e_completion_fail;

  if (data.size() != grid.cells.size())
    RETURN_ERR("bad size %s: %lu %lu", F_DENSITY, data.size(), grid.cells.size());

  for (size_t i = 0; i < data.size(); i++) {
    grid.cells[i].phys_val.d = data[i];
  }
#endif

  return e_completion_success;
}

#endif //! SOLVERS && ILLUM