
#if defined SOLVERS && defined HLLC
#include "hllc_utils.h"

#include "global_value.h"
#include "hllc_flux.h"
#include "linear_alg.h"

#include <omp.h>

Type hllc::GetTimeStep(const hllc_value_t &hllc_set, const std::vector<elem_t> &cells) {

  Type c_max = -1;
  for (auto &el : cells) {
    const Type a = sqrt(kGamma1 * el.phys_val.p / el.phys_val.d);
    const Type c = el.phys_val.v.norm() + a;
    if (c > c_max)
      c_max = c;
  }
  const Type t = hllc_set.CFL * hllc_set.h_min / c_max;
  DIE_IF(t < 0);
  return t;
}

void hllc::HllcPhysToConv(std::vector<elem_t> &cells) {

#pragma omp parallel for
  for (int i = 0; i < cells.size(); i++) {
    GetConvValue(cells[i].phys_val, cells[i].conv_val);
  }
}

void hllc::HllcConvToPhys(std::vector<elem_t> &cells) {

#pragma omp parallel for
  for (int i = 0; i < cells.size(); i++) {
    GetPhysValue(cells[i].conv_val, cells[i].phys_val);
  }
}

void hllc::BoundConditions(const face_t &f, const std::vector<elem_t> &cells, flux_all_t &bound) {
  const int id_l = f.geo.id_l;
  const int id_r = f.geo.id_r;
  const elem_t &cell = cells[id_l];
  Matrix3 T;

  switch (id_r) // id соседа она же признак ГУ
  {
  case e_bound_free:
    bound.conv_val = cell.conv_val;
    break;
  case e_bound_inner_source: {
    bound.phys_val.d = 0.1;
    bound.phys_val.v << 0, 0, 0;
    bound.phys_val.p = 0.1;
    GetConvValue(bound.phys_val, bound.conv_val);
  }

  break;
  case e_bound_out_source:
    bound.phys_val.d = 0.1;
    bound.phys_val.v << 0.2, 0, 0;
    bound.phys_val.p = 0.01;
    GetConvValue(bound.phys_val, bound.conv_val);

    break;
  case e_bound_lock:
    bound.conv_val = cell.conv_val;
    GetRotationMatrix(f.geo.n, T);

    bound.conv_val.v = T * bound.conv_val.v;
    bound.conv_val.v[0] = -bound.conv_val.v[0];
    bound.conv_val.v = (T.transpose()) * bound.conv_val.v;
    break;

  default:
    DIE_IF(f.geo.id_r < 0);

    bound.conv_val = cells[id_r].conv_val;
    break;
  }
}
#endif // HLLC