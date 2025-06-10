#if (defined RHLLC || defined HLLC) && defined USE_MPI
#include "hydro_mpi_cfg.h"

#include "global_value.h"
#include "mpi_shifts.h"
#include "reader_txt.h"

#include <numeric>

using namespace rrhd;

#ifdef RRHD_DEBUG
// #define rhllc_mpi_log(...) WRITE_LOG(__VA_ARGS__)
#define rhllc_mpi_log(...)
#else
#define rhllc_mpi_log(...)
#endif

void hydro_mpi::StartPhysCast(mpi_hd_t *st, grid_t &grid) {
  int np = get_mpi_np(st->comm);
  st->requests_cast_phys.resize(np);

  for (int id = 0; id < np; id++) {
    int N = st->maps[id].cell.right_id - st->maps[id].cell.left_id;
    MPI_Ibcast(grid.cells.data() + st->maps[id].cell.left_id, N,
               MPI_phys_val_t, id, st->comm, &st->requests_cast_phys[id]);
  }
}

void hydro_mpi::SyncPhysCast(mpi_hd_t *st) {
  MPI_Waitall(st->requests_cast_phys.size(), st->requests_cast_phys.data(), MPI_STATUSES_IGNORE);
}

void hydro_mpi::SyncAndCalcPhysCast(grid_t &grid) {
#if 0  
  int myid = get_mpi_id();

  int idx = myid;

  do
  {

    int start = grid.mpi_cfg->disp_cells[idx];
    int end = start + grid.mpi_cfg->send_cells[idx];

#pragma omp parallel default(none) firstprivate(start, end) shared(grid)
    {
#pragma omp for
      for (int i = start; i < end; i++)
      {
        grid.cells[i].cell_data->Init(&grid.cells[i].phys_val);
      }
    }

    MPI_Waitany(grid.mpi_cfg->requests_cast_phys.size(), grid.mpi_cfg->requests_cast_phys.data(), &idx, MPI_STATUSES_IGNORE);
  } while (idx != MPI_UNDEFINED);
#else
  D_LD;
#endif
}

void hydro_mpi::StartExchangeBoundaryCells(mpi_hd_t *st) {
  static bool was_sended = false;
  if (LIKELY(was_sended)) {
    MPI_Waitall(st->requests_send_faces.size(), st->requests_send_faces.data(), MPI_STATUSES_IGNORE);
  }
  MPI_Startall(st->requests_send_faces.size(), st->requests_send_faces.data());
  MPI_Startall(st->requests_rcv_faces.size(), st->requests_rcv_faces.data());
  was_sended = true;
}

void hydro_mpi::SyncExchangeBoundaryCells(mpi_hd_t *st) {
  MPI_Waitall(st->requests_rcv_faces.size(), st->requests_rcv_faces.data(), MPI_STATUSES_IGNORE);
}

void hydro_mpi::InitMpiConfig(const std::vector<int> &metis_id, grid_t &grid, MPI_Comm comm) {

  const int np = get_mpi_np(comm);
  const int myid = get_mpi_id(comm);

  // Вычисление смещений и размера подобластей
  {
    std::vector<int> send;
    std::vector<int> disp;
    SetShifts(metis_id, np, send, disp);

    const int N = grid.size;
    grid.loc_size = send[myid]; // это приведет к правкам на видеокарте(возможно это уже учтено. Надо проверить)
    grid.loc_shift = disp[myid];
  }

  mpi_hd_t *st = new mpi_hd_t;
  st->comm = comm;
  st->maps.resize(np);
  st->maps[myid].cell.left_id = grid.loc_shift;
  st->maps[myid].cell.right_id = grid.loc_shift + grid.loc_size;

  // определение геометрической конфигурации соседей
  // привязка ячейка - узел, грань - узлы
  {
    for (int id_cell = 0; id_cell < grid.size; id_cell++) {
      grid.cells[id_cell].geo.node = metis_id[id_cell];
    }

    for (size_t id_faces = 0; id_faces < grid.faces.size(); id_faces++) {
      face_t *f = &grid.faces[id_faces];
      int idl = f->geo.id_l;
      int idr = f->geo.id_r;

      f->geo.is_regular = 0;
      f->geo.id_l_node = metis_id[idl];

      if (idr < 0) {
        f->geo.id_r_node = idr;
        f->geo.is_regular = (metis_id[idl] == myid);
        continue;
      }

      f->geo.id_r_node = metis_id[idr];

      if ((metis_id[idr] == metis_id[idl]) && (metis_id[idl] == myid)) {
        f->geo.is_regular = 1; // полностью у нас
        continue;
      }
    }
  }

  // Определение соседних ячеек, у которых все грани регулярные
  // т.е. могут быть вычислены из данных на узле
  {

    int left_max = st->maps[myid].cell.left_id;   // самая левая регулярная ячейка
    int right_min = st->maps[myid].cell.right_id; // самая правая регулярная ячейка

    int left_np = -1;  // узел слева
    int right_np = -1; // узел справа

    int cur_left = left_max;
    int cur_right = right_min;

    for (int curCell = cur_left; curCell < cur_right; curCell++) {
      for (int j = 0; j < CELL_SIZE; j++) {
        int f_id = grid.cells[curCell].geo.id_faces[j];
        const geo_face_t *f = &grid.faces[f_id].geo;

        if (!f->is_regular) // эта ячейка не регулярная
        {
          int node = (f->id_l_node != myid) ? f->id_l_node : f->id_r_node;
          if (node >= 0) {
            DIE_IF(abs(node - myid) > 1);

            if (node > myid) {
              if ((right_np == -1) || (right_np == node)) {
                right_min = std::min(curCell, right_min);
                right_np = node;
              }
            } else {
              if ((left_np == -1) || (left_np == node)) {
                left_max = std::max(curCell, left_max);
                left_np = node;
              }
            }
          }
        }
      }
    }

    st->maps[myid].np_l = left_np;
    st->maps[myid].np_r = right_np;

    st->maps[myid].cell.reg_l = left_max + 1; //+1 т.к. left_max - последняя нерегулярная ячейка
    st->maps[myid].cell.reg_r = right_min;    // end()
  }

  // синхронизируем процессы (рассылаем карту всем узлам)
  for (int i = 0; i < np; i++) {
    MPI_Bcast(&st->maps[i], 1, MPI_hllc_map_t, i, comm);
    // MPI_Bcast(st->maps.data() + i, sizeof(st->maps[0]) / sizeof(int64_t), MPI_INT64_T, i, MPI_COMM_WORLD);
  }

  // Определение максимального диапазона граней общитываемого на узле
  {
    int f_min = std::numeric_limits<int>::max();
    int f_max = -1;
    for (size_t id = st->maps[myid].cell.left_id; id < st->maps[myid].cell.right_id; id++) {
      for (size_t j = 0; j < CELL_SIZE; j++) {
        int f_id = grid.cells[id].geo.id_faces[j];
        f_min = std::min(f_min, f_id);
        f_max = std::max(f_max, f_id);
      }
    }

    st->maps[myid].face.left_id = f_min;
    st->maps[myid].face.right_id = f_max + 1; //"за конец"
  }

  // Определение граней, у которых все ячейки регулярные
  // (т.е. могут быть вычислены из данных на узле)
  {
    int f_reg_min = -1;
    int f_reg_max = std::numeric_limits<int>::max();

    int np_l = st->maps[myid].np_l; // узел слева
    int np_r = st->maps[myid].np_r; // узел справа

    // самая левая возможная грань (в нахлест к соседу)
    int left = st->maps[myid].cell.left_id;
    if (np_l >= 0)
      left = st->maps[np_l].cell.reg_r;

    // самая правая возможная грань (в нахлест к соседу)
    int right = st->maps[myid].cell.right_id;
    if (np_r >= 0)
      right = st->maps[np_r].cell.reg_l;

    for (int id = left; id < right; id++) {
      for (int j = 0; j < CELL_SIZE; j++) {
        int f_id = grid.cells[id].geo.id_faces[j];
        const geo_face_t *f = &grid.faces[f_id].geo;

        // соседний узел
        int node = (f->id_l_node != myid) ? f->id_l_node : f->id_r_node;

        if (!f->is_regular && node >= 0) // эта ячейка не регулярная
        {
          if (node == st->maps[myid].np_l) {
            f_reg_min = std::max(f_reg_min, f_id); // подходим слева
          } else {
            f_reg_max = std::min(f_reg_max, f_id); // подходим справа
          }
        }
      }
    }

    // отсекаем по границе глобальной сетки
    f_reg_min = std::max(0, f_reg_min);
    f_reg_max = std::min((int)grid.size_face, f_reg_max);

    st->maps[myid].face.reg_l = f_reg_min + 1; // f_reg_min - последняя не регулярная ячейка
    st->maps[myid].face.reg_r = f_reg_max;
  }

  {
    // синхронизируем процессы (рассылаем карту всем узлам)
    for (int i = 0; i < np; i++) {
      MPI_Bcast(&st->maps[i], 1, MPI_hllc_map_t, i, comm);
      // MPI_Bcast(st->maps.data() + i, sizeof(st->maps[0]) / sizeof(int64_t), MPI_INT64_T, i, MPI_COMM_WORLD);

      rhllc_mpi_log("rcv %d: %d<=>%d, cells: [%d [%d, %d) %d), faces: [%d [%d, %d) %d),\n",
                    i, st->maps[i].np_l, st->maps[i].np_r,
                    st->maps[i].cell.left_id,
                    st->maps[i].cell.reg_l,
                    st->maps[i].cell.reg_r,
                    st->maps[i].cell.right_id,
                    st->maps[i].face.left_id,
                    st->maps[i].face.reg_l,
                    st->maps[i].face.reg_r,
                    st->maps[i].face.right_id);
    }
  }

  // инициализируем MPI запросы
  {
    st->requests_send_faces.clear();
    st->requests_rcv_faces.clear();
    st->requests_cast_phys.clear();

    int np_lr[2] = {static_cast<int>(st->maps[myid].np_l),
                    static_cast<int>(st->maps[myid].np_r)};

    // Ячейка лежат: { [left_id [reg_l, reg_r), right_id) }
    int lr_min[2] = {static_cast<int>(st->maps[myid].cell.left_id),
                     static_cast<int>(st->maps[myid].cell.reg_r)};

    int lr_len[2] = {static_cast<int>(st->maps[myid].cell.reg_l - lr_min[0]),
                     static_cast<int>(st->maps[myid].cell.right_id - lr_min[1])};

    for (int i = 0; i < 2; i++) {
      if (lr_len[i] > 0 && np_lr[i] >= 0) {
        MPI_Request rq;
        MPI_Send_init(&grid.cells[lr_min[i]], lr_len[i], MPI_flux_elem_t, np_lr[i], i,
                      st->comm, &rq);
        st->requests_send_faces.push_back(rq);

        rhllc_mpi_log("send: N=%d, adr=%d, src=%d dst=%d, tag=%d\n",
                      lr_len[i], lr_min[i], myid, np_lr[i], i);
      }
    }

    // сосед слева
    if (np_lr[0] >= 0) {
      int l_min = st->maps[np_lr[0]].cell.reg_r;
      int l_len = st->maps[np_lr[0]].cell.right_id - l_min;

      if (l_len > 0) {
        MPI_Request rq;
        MPI_Recv_init(&grid.cells[l_min], l_len, MPI_flux_elem_t, np_lr[0], 1,
                      st->comm, &rq);
        st->requests_rcv_faces.push_back(rq);
        rhllc_mpi_log("rcv: N=%d, adr=%d, src=%d, dst=%d, tag=%d\n",
                      l_len, l_min, np_lr[0], myid, 1);
      }
    }

    // cосед справа
    if (np_lr[1] >= 0) {
      int r_min = st->maps[np_lr[1]].cell.left_id;
      int r_len = st->maps[np_lr[1]].cell.reg_l - r_min;
      if (r_len > 0) {
        MPI_Request rq;
        MPI_Recv_init(&grid.cells[r_min], r_len, MPI_flux_elem_t, np_lr[1], 0,
                      st->comm, &rq);
        st->requests_rcv_faces.push_back(rq);
        rhllc_mpi_log("rcv: N=%d, adr=%d, src=%d,dst=%d, tag=%d\n",
                      r_len, r_min, np_lr[1], myid, 0);
      }
    }
  }

  if (grid.mpi_cfg) {
    delete grid.mpi_cfg;
  }
  grid.mpi_cfg = st;
}
