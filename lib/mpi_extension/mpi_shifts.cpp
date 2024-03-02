#include "mpi_shifts.h"

void GetSend(const int np, const IdType n, std::vector<IdType> &send_count) {

  send_count.assign(np, n / np);

  if (n % np) // если число процессов не кратно размерности задачи
  {
    for (int i = 0; i < n % np; i++) // первые процессы берут на единицу больше тел
      ++send_count[i];
  }
}

void GetDisp(const int np, const IdType n, std::vector<IdType> &disp) {

  disp.assign(np, 0);
  int a = n % np;
  for (int i = 1; i < np; i++)
    disp[i] = i * (n / np) + a;

  if (n % np) // если число процессов не кратно размерности задачи
  {
    for (int i = 1; i < a; i++) // смещения для процессов за ними увеличивается начиная со второго
      disp[a - i] -= i;
  }
}
#ifdef USE_MPI
#include "dbgdef.h"
#include <algorithm>

void SetShifts(const std::vector<int> &metis_id, mpi_hllc_t *cfg) {
  int np = get_mpi_np(cfg->comm);

  int np_count = *std::max_element(metis_id.begin(), metis_id.end());
  DIE_IF(((np_count + 1) != np));

  cfg->send_cells.resize(np, 0);
  cfg->disp_cells.resize(np, 0);
  for (int i = 0; i < np; i++) {
    cfg->send_cells[i] = std::count(metis_id.begin(), metis_id.end(), i);

    cfg->disp_cells[i] = 0;
    for (int j = 0; j < i; j++) {
      cfg->disp_cells[i] += cfg->send_cells[j];
    }
  }
}
#endif