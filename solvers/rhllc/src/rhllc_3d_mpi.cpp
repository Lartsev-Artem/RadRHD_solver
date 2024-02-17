#if defined RHLLC && defined SOLVERS && USE_MPI
#include "rhllc_mpi.h"

#include "global_value.h"
#include "reader_txt.h"

#include "mpi_ext.h"
#include "mpi_shifts.h"

#include "rhllc_flux_stab.h"
#include "rhllc_utils.h"

using namespace rhllc;

struct mpi_hllc_t {
  std::vector<IdType> send_cells;        ///< кол-во отправок для каждого процесса
  std::vector<IdType> disp_cells;        ///< смещения по ячейкам для процессов
  std::vector<IntId> id_irregular_faces; ///< номера граней с границе областей процессов

  std::vector<MPI_Request> requests_cast_phys;  ///< запросы на передачу физ. переменных
  std::vector<MPI_Request> requests_send_faces; ///< запросы на передачу потоков грани
  std::vector<MPI_Request> requests_rcv_faces;  ///< запросы на приём потоков грани
};

void InitCast(mpi_hllc_t &hllc_st, grid_t &grid, const MPI_Comm &comm) {

  int np = get_mpi_np(comm);
  hllc_st.requests_cast_phys.resize(np);

  for (int id = 0; id < np; id++) {
    MPI_Ibcast(grid.cells.data() + hllc_st.disp_cells[id], hllc_st.send_cells[id], MPI_phys_val_t, id, comm, &hllc_st.requests_cast_phys[id]);
  }
}

void SyncCast(const mpi_hllc_t &hllc_st) {
  MPI_Waitall(hllc_st.requests_cast_phys.size(), hllc_st.requests_cast_phys.data(), MPI_STATUSES_IGNORE);
  // MPI_Waitany( int count , MPI_Request array_of_requests[] , int* index , MPI_Status* status);
}

void start_flux(const mpi_hllc_t &hllc_st) {
  static bool was_sended = false;
  if (LIKELY(was_sended)) {
    MPI_Waitall(hllc_st.requests_send_faces.size(), hllc_st.requests_send_faces.data());
  }
  MPI_Startall(hllc_st.requests_send_faces.size(), hllc_st.requests_send_faces.data());
  MPI_Startall(hllc_st.requests_rcv_faces.size(), hllc_st.requests_rcv_faces.data());
  was_sended = true;
}

void wait_flux(const mpi_hllc_t &hllc_st) {
  MPI_Waitall(hllc_st.requests_rcv_faces.size(), hllc_st.requests_rcv_faces.data());
}

void Init(mpi_hllc_t &hllc_st, grid_t &grid, const MPI_Comm &comm) {

  int np = get_mpi_np(comm);
  int myid = get_mpi_id(comm);

  const int N = grid.size;

  GetSend(np, N, hllc_st.send_cells);
  GetDisp(np, N, hllc_st.disp_cells);

  grid.loc_size = hllc_st.send_cells[myid]; //это приведет к правкам на видеокарте(возможно это уже учтено. Надо проверить)

  std::vector<int> metis_id;
  files_sys::txt::ReadData(glb_files.base_address + F_SEPARATE_METIS(np), metis_id);

  hllc_st.id_irregular_faces.clear();
  hllc_st.requests_send_faces.clear();
  hllc_st.requests_rcv_faces.clear();
  hllc_st.requests_cast_phys.clear();

  for (int id_cell = 0; id_cell < grid.size; id_cell++) {
    grid.cells[id_cell].geo.node = metis_id[id_cell];
  }

  hllc_st.id_irregular_faces.reserve(N / np);

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
      f->geo.is_regular = 1; //полностью у нас
      continue;
    }

    if (metis_id[idl] == myid) {
      MPI_Request rq_send;
      MPI_Request rq_rcv;

      MPI_Send_init(&grid.cells[idl], 1, MPI_flux_elem_t, metis_id[idr], 0, comm, &rq_send);
      MPI_Recv_init(&grid.cell[idl], 1, MPI_flux_elem_t, metis_id[idr], 0, comm, &rq_rcv);

      hllc_st.requests_send_faces.push_back(rq_send);
      hllc_st.requests_rcv_faces.push_back(rq_rcv);
      hllc_st.id_irregular_faces.push_back(id_faces); // только на нашем узле
    }

    if (metis_id[idr] == myid) {
      MPI_Request rq_send;
      MPI_Request rq_rcv;

      MPI_Send_init(&grid.cells[idr], 1, MPI_flux_elem_t, metis_id[idl], 0, comm, &rq_send);
      MPI_Recv_init(&grid.cell[idr], 1, MPI_flux_elem_t, metis_id[idl], 0, comm, &rq_rcv);

      requests_send_faces.push_back(rq_send);
      requests_rcv_faces.push_back(rq_rcv);
      id_irregular_faces.push_back(id_faces); // только на нашем узле
    }
  }
  id_irregular_faces.shrink_to_fit();
}

void rhllc_mpi::Hllc3dStab(const Type tau, grid_t &grid) {
  rhllc::max_signal_speed = 0;

#pragma omp parallel default(none) firstprivate(tau) shared(grid, glb_files, rhllc::max_signal_speed)
  {
    const int size_grid = grid.size;
    Type max_speed = 0;
    flux_all_t bound_val;
    const int size_face = grid.faces.size();
// потоки
#pragma omp for
    for (int i = 0; i < size_face; i++) {
      face_t &f = grid.faces[i];
      BoundConditions(f, grid.cells, bound_val);
      max_speed = std::max(max_speed, GetFluxStab(grid.cells[f.geo.id_l].conv_val, bound_val.conv_val, grid.cells[f.geo.id_l].phys_val, bound_val.phys_val, f));
    }

#pragma omp for
    for (int i = 0; i < size_grid; i++) {
      elem_t &el = grid.cells[i];
      flux_t sumF;
      for (int j = 0; j < CELL_SIZE; j++) {
        if (el.geo.sign_n[j]) {
          sumF += grid.faces[el.geo.id_faces[j]].f;
        } else {
          sumF -= grid.faces[el.geo.id_faces[j]].f;
        }
      }
      sumF *= (tau / el.geo.V);
      el.conv_val -= sumF;

      if (GetPhysValueStab(el.conv_val, el.phys_val)) {
        DIE_IF(PhysPressureFix(el.conv_val, el.phys_val));
      }
    }

    if (max_speed > rhllc::max_signal_speed) {
#pragma omp critical
      {
        rhllc::max_signal_speed = std::max(max_signal_speed, max_speed);
      }
    }

  } // omp
}

#endif

#if 1

void hllc() {
  int N = 9;
  int np = 3;
  int metis_id[9] = {0, 0, 0, 1, 1, 1, 2, 2, 2};
  int id_irreg_faces[100];

  // (if != 0) wait_send_all(); //завершаем предыдущие вызовы
  // start_all sendrcv

  for (auto &f : face) {
    if (f.is_reg) // (f.idr==f.idl && f_idr==myid)
    {
      // calc
    }
  }

  // wait_rcv_all

  for (auto id : id_irreg_faces) {
    // calc: faces[id]
  }

  for (auto c : cells) {
    if (c.node == myid) {
      // calc
    }
  }

  for (int id = 0; id < np; id++) {
    IBcast(cells.data() + disp[id], send[id], id, MPI_COMM_WORLD, &requests_cast_phys[id]);
  }

  //здесь можно сделать  так:

  int rcvCast = np - 1; //кроме себя
  do {
    flag = false;
    for (int id = 0; id < np; id++) {
      if (myid != id && requests_cast_phys[id] != NULL)
        if (Mpi_test(requests_cast_phys[id])) {
          // calc(alpha,betta,T,logT)[disp[id]:send[id]+disp[id]]
          --rcvCast;
          requests_cast_phys[id] != NULL;
        }
    }
  } while (rcvCast);

  Mpiwaitall(requests_cast_phys.size(), requests_cast_phys.data());
}
#endif