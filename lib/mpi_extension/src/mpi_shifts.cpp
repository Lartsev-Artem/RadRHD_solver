#include "mpi_shifts.h"

void GetSend(const int np, const IdType n, std::vector<IdType> &send_count)
{

  send_count.assign(np, n / np);

  if (n % np) // если число процессов не кратно размерности задачи
  {
    for (int i = 0; i < n % np; i++) // первые процессы берут на единицу больше тел
      ++send_count[i];
  }
}

void GetDisp(const int np, const IdType n, std::vector<IdType> &disp)
{

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

void SetShifts(const std::vector<int> &metis_id, int np,
               std::vector<int> &send, std::vector<int> &disp)
{
  int np_count = *std::max_element(metis_id.begin(), metis_id.end());
  DIE_IF(((np_count + 1) != np));

  send.resize(np, 0);
  disp.resize(np, 0);
  for (int i = 0; i < np; i++)
  {
    send[i] = std::count(metis_id.begin(), metis_id.end(), i);

    disp[i] = 0;
    for (int j = 0; j < i; j++)
    {
      disp[i] += send[j];
    }
  }
}
#endif