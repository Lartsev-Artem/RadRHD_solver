/**
 * @file foo.cpp
 * @brief Файл ни куда не подключается. Присутствует временно
   @todo радиационное ускорение
 */

#include "solvers_struct.h"

///\todo при переходе к релятивистской постановке пропали коэффициенты 1/c и 1/c^2
int SolveIllumAndHLLC(const Type tau, grid_t &grid) {
  Type min = 1;
  int ret_flag = 0;

#pragma omp parallel default(none) firstprivate(tau) shared(ret_flag, grid)
  {
    const int n = grid.size;
#pragma omp for
    for (int i = 0; i < n; i++) {
      elem_t &el = grid.cells[i];

#ifdef USE_CUDA
      el.conv_val.v += tau * (-grid.divimpuls[i] * 1); //- tau*(stream[i][0] - prev_stream[i][0]) / tau / c / c);  // +F_g //vx
      el.conv_val.p += tau * (-grid.divstream[i]);     //  tau*(-(energy[i] - prev_energy[i]) / tau / c)  //+F_g.dot(vel)  //e

#else
      el.conv_val.v += tau * (-el.illum_val.div_impuls * 1); //- tau*(stream[i][0] - prev_stream[i][0]) / tau / c / c);  // +F_g //vx
      el.conv_val.p += tau * (-el.illum_val.div_stream);     //  tau*(-(energy[i] - prev_energy[i]) / tau / c)  //+F_g.dot(vel)  //e
#endif
#ifdef DEBUG
      if (el.conv_val.p - sqrt(el.conv_val.d * el.conv_val.d + el.conv_val.v.dot(el.conv_val.v)) < 0) // условие физичности консервативных переменных
      {
        ret_flag = 1; // return 1;
      }
#endif
    }
  }
  return ret_flag;
}
