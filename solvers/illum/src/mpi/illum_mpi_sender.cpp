#if defined ILLUM && defined SOLVERS && defined USE_MPI
#include "illum_mpi_sender.h"
#include "illum_mpi_struct.h"
#include "mpi_shifts.h"

MPI_Comm MPI_COMM_ILLUM = MPI_COMM_WORLD;
illum::mpi_sender_t section_1;

#if defined SEPARATE_GPU && !defined SPECTRUM
void illum::separate_gpu::InitSender(const MPI_Comm &comm, const grid_directions_t &grid_dir, const grid_t &grid) {

  MpiInitStruct(grid_dir);
  int np = get_mpi_np(comm);
  int myid = get_mpi_id(comm);

  std::vector<IdType> disp;
  std::vector<IdType> send_count;

  GetSend(np, grid_dir.size, send_count);
  GetDisp(np, grid_dir.size, disp);

  const IdType size_section = send_count[myid];
  section_1.size = size_section;
  {
    const IdType N = grid_dir.size - size_section;
    section_1.requests_rcv.resize(N, MPI_REQUEST_NULL);
    section_1.status_rcv.resize(N);
    section_1.flags_send_to_gpu.resize(N, 0);

    IdType size_msg = grid.size;
    DIE_IF(size_msg > ((1u << 31) - 1));

    int cc = 0;
    for (int src = 0; src < np; src++) {
      if (src == myid)
        continue;
      for (int j = 0; j < send_count[src]; j++) {
        int tag = disp[src] + j;
        MPI_Recv_init(grid.Illum + tag, (int)size_msg, MPI_RECV_ILLUM_T, src, tag, comm, &section_1.requests_rcv[cc++]);
      }
    }

    section_1.requests_send.resize(size_section * (np - 1), MPI_REQUEST_NULL);

    for (int num_direction = 0; num_direction < size_section; num_direction++) {
      const int tag = disp[myid] + num_direction; // teg соответствует номеру направления
      cc = 0;
      for (int id = 0; id < np; id++) {
        if (id == myid)
          continue;

        MPI_Send_init(grid.Illum + tag, (int)size_msg, MPI_RECV_ILLUM_T, id, tag, comm, &section_1.requests_send[(np - 1) * num_direction + cc++]);
      }
    }
  }

  WRITE_LOG("Init mpi sender\n");
  return;
}
#endif

#ifdef SPECTRUM
void illum::spectrum_gpu::InitSender(const MPI_Comm &comm, const grid_directions_t &grid_dir, const grid_t &grid) {

  // MPI_Send(Illum_local.data(), N*F,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
  // MPI_Recv(Illum.data(), N,subarray2,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  spectrum::MpiInitStruct(grid);

  int np = get_mpi_np(comm);
  int myid = get_mpi_id(comm);

  std::vector<IdType> disp;
  std::vector<IdType> send_count;

  GetSend(np, grid_dir.size, send_count);
  GetDisp(np, grid_dir.size, disp);

  const IdType size_section = send_count[myid];
  section_1.size = size_section;
  {
    const IdType N = grid_dir.size - size_section;
    section_1.requests_rcv.resize(N, MPI_REQUEST_NULL);
    section_1.status_rcv.resize(N);
    section_1.flags_send_to_gpu.resize(N, 0);

    IdType size_msg = grid.size;
    DIE_IF(size_msg > ((1u << 31) - 1));

    int cc = 0;
    for (int src = 0; src < np; src++) {
      if (src == myid)
        continue;
      for (int j = 0; j < send_count[src]; j++) {
        int tag = disp[src] + j;
        MPI_Recv_init(grid.Illum + tag * grid.size_frq, (int)size_msg, MPI_SUB_ARRAY_RECV_T, src, tag, comm, &section_1.requests_rcv[cc++]);
      }
    }

    section_1.requests_send.resize(size_section * (np - 1), MPI_REQUEST_NULL);

    for (int num_direction = 0; num_direction < size_section; num_direction++) {
      const int tag = disp[myid] + num_direction; // teg соответствует номеру направления
      cc = 0;
      for (int id = 0; id < np; id++) {
        if (id == myid)
          continue;

        MPI_Send_init(grid.Illum + tag * grid.size_frq, (int)size_msg, MPI_SUB_ARRAY_RECV_T, id, tag, comm, &section_1.requests_send[(np - 1) * num_direction + cc++]);
      }
    }
  }

  WRITE_LOG("Init mpi sender\n");
  return;
}
#endif

illum::mpi_sender_t section_2;
std::vector<IdType> disp_illum;
void illum::gpu_async::InitSender(const grid_directions_t &grid_dir, const grid_t &grid) {

  int np = get_mpi_np();
  int myid = get_mpi_id();

  std::vector<IdType> disp;
  std::vector<IdType> send_count;

  GetSend(np, grid_dir.size, send_count);
  GetDisp(np, grid_dir.size, disp);

  disp_illum.resize(np);
  for (int i = 0; i < np; i++) {
    disp_illum[i] = disp[i] * CELL_SIZE * grid.size;
  }

  const IdType size_first_section = std::max(((1u << 31) / (CELL_SIZE * grid.size)) - 1, 1 * send_count[0] / 3); // первый узел будем потенциально разгружать(на всех узла должно хватать направлений)
  section_1.size = size_first_section;
  {
    //==================first rcv ==============================
    section_1.requests_rcv.resize(np - 1, MPI_REQUEST_NULL);
    section_1.status_rcv.resize(np - 1);
    section_1.flags_send_to_gpu.resize(np - 1, 0);
    const IdType size_msg = grid.size * CELL_SIZE * size_first_section;

    DIE_IF(size_msg > ((1u << 31) - 1));

    int cc = 0;
    for (int src = 0; src < np; src++) {
      if (src == myid)
        continue;
      int tag = src;
      MPI_Recv_init(grid.Illum + disp_illum[tag] /*size_msg * tag*/, (int)size_msg, MPI_DOUBLE, src, tag, MPI_COMM_ILLUM, &section_1.requests_rcv[cc++]);
    }

    //==================first send ==============================

    cc = 0;
    section_1.requests_send.resize(np - 1, MPI_REQUEST_NULL);
    for (int id = 0; id < np; id++) {
      if (id == myid)
        continue;
      MPI_Send_init(grid.Illum + disp_illum[myid] /*size_msg * myid*/, (int)size_msg, MPI_DOUBLE, id, myid, MPI_COMM_ILLUM, &section_1.requests_send[cc++]);
    }
  }

  //======================MPI_INIT=========================

  const IdType local_size = send_count[myid];
  const IdType size_second_section = local_size - size_first_section;
  section_2.size = size_second_section;
  {
    const IdType N = grid_dir.size - (size_first_section * np) - size_second_section;
    section_2.requests_rcv.resize(N, MPI_REQUEST_NULL);
    section_2.status_rcv.resize(N);
    section_2.flags_send_to_gpu.resize(N, 0);

    IdType size_msg = grid.size * CELL_SIZE;
    DIE_IF(size_msg > ((1u << 31) - 1));

    int cc = 0;
    for (int src = 0; src < np; src++) {
      if (src == myid)
        continue;
      for (int j = size_first_section; j < send_count[src]; j++) {
        int tag = disp[src] + j;
        MPI_Recv_init(grid.Illum + size_msg * tag, (int)size_msg, MPI_DOUBLE, src, tag, MPI_COMM_ILLUM, &section_2.requests_rcv[cc++ /*tag - local_size*/]);
      }
    }

    section_2.requests_send.resize(size_second_section * (np - 1), MPI_REQUEST_NULL);

    for (int num_direction = size_first_section; num_direction < local_size; num_direction++) {
      const int tag = disp[myid] + num_direction; // teg соответствует номеру направления
      cc = 0;
      for (int id = 0; id < np; id++) {
        if (id == myid)
          continue;

        MPI_Send_init(grid.Illum + (disp[myid] + num_direction) * size_msg, (int)size_msg, MPI_DOUBLE, id, tag, MPI_COMM_ILLUM, &section_2.requests_send[(np - 1) * (num_direction - size_first_section) + cc++]);
      }
    }
  }

  WRITE_LOG("Init mpi sender\n");
  return;
}

#endif //! defined ILLUM && defined SOLVERS  && !defined USE_MPI