#if defined RHLLC && defined SOLVERS && defined USE_MPI
#include "rhllc_mpi.h"

#include "global_value.h"
#include "reader_txt.h"

#include "mpi_shifts.h"

#include "rhllc_flux_stab.h"
#include "rhllc_utils.h"

void rhllc_mpi::StartPhysCast(mpi_hllc_t &hllc_st, grid_t &grid) {

  int np = get_mpi_np(hllc_st.comm);
  hllc_st.requests_cast_phys.resize(np);

  for (int id = 0; id < np; id++) {
    MPI_Ibcast(grid.cells.data() + hllc_st.disp_cells[id], hllc_st.send_cells[id],
               MPI_phys_val_t, id, hllc_st.comm, &hllc_st.requests_cast_phys[id]);
  }
}

void rhllc_mpi::SyncPhysCast(mpi_hllc_t &hllc_st) {
  MPI_Waitall(hllc_st.requests_cast_phys.size(), hllc_st.requests_cast_phys.data(), MPI_STATUSES_IGNORE);
}

void rhllc_mpi::SyncAndCalcPhysCast(grid_t &grid) {
  int myid = get_mpi_id();

  int idx = myid;

  do {

    int start = grid.mpi_cfg->disp_cells[idx];
    int end = start + grid.mpi_cfg->send_cells[idx];

#pragma omp parallel default(none) firstprivate(start, end) shared(grid)
    {
#pragma omp for
      for (int i = start; i < end; i++) {
        grid.cells[i].cell_data->Init(&grid.cells[i].phys_val);
      }
    }

    MPI_Waitany(grid.mpi_cfg->requests_cast_phys.size(), grid.mpi_cfg->requests_cast_phys.data(), &idx, MPI_STATUSES_IGNORE);
  } while (idx != MPI_UNDEFINED);
}

void rhllc_mpi::StartExchangeBoundaryCells(mpi_hllc_t &hllc_st) {
  static bool was_sended = false;
  if (LIKELY(was_sended)) {
    MPI_Waitall(hllc_st.requests_send_faces.size(), hllc_st.requests_send_faces.data(), MPI_STATUSES_IGNORE);
  }
  MPI_Startall(hllc_st.requests_send_faces.size(), hllc_st.requests_send_faces.data());
  MPI_Startall(hllc_st.requests_rcv_faces.size(), hllc_st.requests_rcv_faces.data());
  was_sended = true;
}

void rhllc_mpi::SyncExchangeBoundaryCells(mpi_hllc_t &hllc_st) {
  MPI_Waitall(hllc_st.requests_rcv_faces.size(), hllc_st.requests_rcv_faces.data(), MPI_STATUSES_IGNORE);
}

void rhllc_mpi::InitMpiConfig(const std::vector<int> &metis_id, grid_t &grid, mpi_hllc_t *hllc_st) {

  hllc_st->comm = MPI_COMM_WORLD;

  int np = get_mpi_np(hllc_st->comm);
  int myid = get_mpi_id(hllc_st->comm);

  const int N = grid.size;

  GetSend(np, N, hllc_st->send_cells);
  GetDisp(np, N, hllc_st->disp_cells);

  grid.loc_size = hllc_st->send_cells[myid]; //это приведет к правкам на видеокарте(возможно это уже учтено. Надо проверить)

  hllc_st->id_irregular_faces.clear();
  hllc_st->requests_send_faces.clear();
  hllc_st->requests_rcv_faces.clear();
  hllc_st->requests_cast_phys.clear();

  for (int id_cell = 0; id_cell < grid.size; id_cell++) {
    grid.cells[id_cell].geo.node = metis_id[id_cell];
  }

  hllc_st->id_irregular_faces.reserve(N / np);

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

      MPI_Send_init(&grid.cells[idl], 1, MPI_flux_elem_t, metis_id[idr], 0, hllc_st->comm, &rq_send);
      MPI_Recv_init(&grid.cells[idl], 1, MPI_flux_elem_t, metis_id[idr], 0, hllc_st->comm, &rq_rcv);

      hllc_st->requests_send_faces.push_back(rq_send);
      hllc_st->requests_rcv_faces.push_back(rq_rcv);
      hllc_st->id_irregular_faces.push_back(id_faces); // только на нашем узле
    }

    if (metis_id[idr] == myid) {
      MPI_Request rq_send;
      MPI_Request rq_rcv;

      MPI_Send_init(&grid.cells[idr], 1, MPI_flux_elem_t, metis_id[idl], 0, hllc_st->comm, &rq_send);
      MPI_Recv_init(&grid.cells[idr], 1, MPI_flux_elem_t, metis_id[idl], 0, hllc_st->comm, &rq_rcv);

      hllc_st->requests_send_faces.push_back(rq_send);
      hllc_st->requests_rcv_faces.push_back(rq_rcv);
      hllc_st->id_irregular_faces.push_back(id_faces); // только на нашем узле
    }
  }
  hllc_st->id_irregular_faces.shrink_to_fit();
}

void rhllc_mpi::Hllc3dStab(const Type tau, grid_t &grid) {

  rhllc::max_signal_speed = 0;
  MPI_Request rq_max_speed = MPI_REQUEST_NULL;
  int myid = get_mpi_id();

  StartExchangeBoundaryCells(*grid.mpi_cfg);

#pragma omp parallel default(none) firstprivate(tau, myid) shared(grid, glb_files, rhllc::max_signal_speed, rq_max_speed)
  {
    const int size_grid = grid.size;
    Type max_speed = 0;
    flux_all_t bound_val;
    const int size_face = grid.faces.size();
// потоки
#pragma omp for
    for (int i = 0; i < size_face; i++) {
      face_t &f = grid.faces[i];
      if (f.geo.is_regular) {
        rhllc::BoundConditions(f, grid.cells, bound_val);
        max_speed = std::max(max_speed, rhllc::GetFluxStab(grid.cells[f.geo.id_l].conv_val, bound_val.conv_val, grid.cells[f.geo.id_l].phys_val, bound_val.phys_val, f));
      }
    }

#pragma omp master
    {
      SyncExchangeBoundaryCells(*grid.mpi_cfg);
    }

// потоки на граничных ячейках узла
#pragma omp for
    for (int i = 0; i < grid.mpi_cfg->id_irregular_faces.size(); i++) {
      face_t &f = grid.faces[i];
      rhllc::BoundConditions(f, grid.cells, bound_val);
      max_speed = std::max(max_speed, rhllc::GetFluxStab(grid.cells[f.geo.id_l].conv_val, bound_val.conv_val, grid.cells[f.geo.id_l].phys_val, bound_val.phys_val, f));
    }

    // ищем максимум по потокам
    if (max_speed > rhllc::max_signal_speed) {
#pragma omp critical
      {
        rhllc::max_signal_speed = std::max(rhllc::max_signal_speed, max_speed);
      }
    }
  }
  // ищем максимум по всем узлам
  //#pragma omp master
  {
    Type max_speed = rhllc::max_signal_speed;
    MPI_Iallreduce(&max_speed, &rhllc::max_signal_speed, 1, MPI_DOUBLE, MPI_MAX, grid.mpi_cfg->comm, &rq_max_speed);
  }

#pragma omp parallel default(none) firstprivate(tau, myid) shared(grid, glb_files, rhllc::max_signal_speed, rq_max_speed)
  {
    int start = grid.mpi_cfg->disp_cells[myid];
    int end = start + grid.mpi_cfg->send_cells[myid];

#pragma omp for
    for (int i = start; i < end; i++) {
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

      if (rhllc::GetPhysValueStab(el.conv_val, el.phys_val)) {
        DIE_IF(rhllc::PhysPressureFix(el.conv_val, el.phys_val));
      }
    }

  } // omp

  MPI_Wait(&rq_max_speed, MPI_STATUS_IGNORE);
}

void rhllc_mpi::HllcConvToPhys(grid_t &grid) {
  int myid = get_mpi_id();
  int start = grid.mpi_cfg->disp_cells[myid];
  int end = start + grid.mpi_cfg->send_cells[myid];

#pragma omp parallel default(none) firstprivate(start, end) shared(grid)
  {
#pragma omp for
    for (int i = start; i < end; i++) {
      rhllc::GetPhysValueStab(grid.cells[i].phys_val, grid.cells[i].conv_val);
    }
  }
}

#include "radRHD_utils.h"
void rhllc_mpi::AddRadFlux(grid_t &grid) {

  int myid = get_mpi_id();
  int start = grid.mpi_cfg->disp_cells[myid];
  int end = start + grid.mpi_cfg->send_cells[myid];

  constexpr Type ds = 1.0;
  Type tau = ds * _hllc_cfg.tau;
#pragma omp parallel default(none) firstprivate(start, end, tau) shared(grid)
  {

#pragma omp for
    for (int cell = start; cell < end; cell++) {

      Vector4 G;
      rad_rhd::GetRadSourceOpt(cell, grid, G);
      G *= tau;

      flux_t &conv_val = grid.cells[cell].conv_val;
      conv_val.p += G[0];
      conv_val.v[0] += G[1];
      conv_val.v[1] += G[2];
      conv_val.v[2] += G[3];
    }
  }
}

#endif //!  RHLLC &&  SOLVERS && USE_MPI