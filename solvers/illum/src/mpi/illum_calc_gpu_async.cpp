#if defined ILLUM && defined SOLVERS && defined USE_MPI && defined USE_CUDA
#include "illum_calc_gpu_async.h"

#include "illum_params.h"
#include "illum_utils.h"
#include "scattering.h"

#include "cuda_interface.h"
#include "solvers_struct.h"

#include "mpi_ext.h"
#include "mpi_shifts.h"

#include <omp.h>

#include <chrono>
namespace tick = std::chrono;

#include "illum_mpi_sender.h"
#if 0
struct mpi_sender_t {
  IdType size;                            ///< размер секции по направлениям (не по запросам)
  std::vector<MPI_Request> requests_rcv;  //все запросы сообщений отправки и принятия
  std::vector<MPI_Status> status_rcv;     //статусы всех обменов
  std::vector<int> flags_send_to_gpu;     //флаги указывающие на отправку пакета на gpu
  std::vector<MPI_Request> requests_send; //все запросы сообщений отправки и принятия
};

static std::vector<IdType> disp_illum;
static mpi_sender_t section_1;
static mpi_sender_t section_2;
static MPI_Comm MPI_COMM_ILLUM = MPI_COMM_WORLD;

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

  return;
}
#endif
int illum::gpu_async::CalculateIllum(const grid_directions_t &grid_direction, const std::vector<std::vector<bits_flag_t>> &face_states,
                                     const std::vector<IntId> &neighbours, const std::vector<std::vector<IntId>> &inner_bound_code,
                                     const std::vector<std::vector<cell_local>> &vec_x0, std::vector<BasePointTetra> &vec_x,
                                     const std::vector<std::vector<IntId>> &sorted_id_cell, grid_t &grid) {

  const IdType n_illum = grid.size * CELL_SIZE;

  int np = get_mpi_np();
  int myid = get_mpi_id();

  DIE_IF(np <= 1); //на одном узле не работает. Тогда надо закрывать пересылки mpi

  const IdType local_size = grid_direction.loc_size;
  const IdType local_disp = grid_direction.loc_shift;

  int iter = 0;    ///< номер итерации
  double norm = 0; ///< норма ошибки
  std::vector<Type> norms(np, -10);

  do {
    auto start_clock = tick::steady_clock::now();
    norm = -1;

    if (section_1.requests_send[0] != MPI_REQUEST_NULL) {
      MPI_Waitall(section_1.requests_send.size(), section_1.requests_send.data(), MPI_STATUSES_IGNORE);
    }
    if (section_2.requests_send[0] != MPI_REQUEST_NULL) {
      MPI_Waitall(section_2.requests_send.size(), section_2.requests_send.data(), MPI_STATUSES_IGNORE);
    }
    cuda::interface::CudaSyncStream(cuda::e_cuda_scattering_1);

    if (myid == 0) {
      //самый высоких приоритет, т.к. надо расчитать, до конфликта с асинхронной отправкой
      cuda::interface::CalculateAllParamAsync(grid_direction, grid, cuda::e_cuda_params); //запустим расчёт параметров здесь
                                                                                          // на выходе получим ответ за 1 шаг до сходимости, но зато без ожидания на выходе
    }

    /*---------------------------------- далее FOR по направлениям----------------------------------*/
    const IdType count_directions = grid_direction.size;

#pragma omp parallel default(none) firstprivate(count_directions, n_illum, myid, np, local_disp, local_size)     \
    shared(sorted_id_cell, neighbours, face_states, vec_x0, vec_x, grid, norm, disp_illum, section_1, section_2, \
           inner_bound_code)
    {
#pragma omp single
      {
        section_1.flags_send_to_gpu.assign(section_1.flags_send_to_gpu.size(), 0);
        if (np > 1) {
          MPI_Startall(section_1.requests_rcv.size(), section_1.requests_rcv.data());
        }
      }

      const int count_th = omp_get_num_threads();
      const int num_th = omp_get_thread_num();

      const int count_cells = grid.size;

      Type loc_norm = -1;
      std::vector<Vector3> *inter_coef = &grid.inter_coef_all[omp_get_thread_num()]; ///< указатель на коэффициенты интерполяции по локальному для потока направлению

#pragma omp for
      for (IdType num_direction = 0; num_direction < section_1.size; ++num_direction) {

        const cell_local *X0_ptr = vec_x0[num_direction].data(); ///< индексация по массиву определяющих гранях (конвеерная т.к. заранее не известны позиции точек)
        const IntId *code_bound = inner_bound_code[num_direction].data();
        /*---------------------------------- далее FOR по ячейкам----------------------------------*/
        for (IdType h = 0; h < count_cells; ++h) {

          const IdType num_cell = sorted_id_cell[num_direction][h];
          const IdType face_block_id = num_cell * CELL_SIZE;

          elem_t *cell = &grid.cells[num_cell];

          // расчитываем излучения на выходящих гранях
          for (ShortId num_out_face = 0; num_out_face < CELL_SIZE; ++num_out_face) {

            const IdType neigh_id = neighbours[face_block_id + num_out_face]; ///< сосед к текущей грани

            // если эта грань входящая и граничная, то пропускаем её
            if (CHECK_BIT(face_states[num_direction][num_cell], num_out_face) == e_face_type_in) {
              if (neigh_id < 0) {
#ifdef USE_TRACE_THROUGH_INNER_BOUNDARY
                Type I0;
                if (neigh_id == e_bound_inner_source) {
                  Vector3 I_def = *code_bound >= 0 ? (*inter_coef)[*code_bound] : Vector3::Zero(); //т.е. определять будет ячейка, а не геометрия
                  I0 = illum::BoundaryConditions(neigh_id, *code_bound, I_def);
                  code_bound++;
                } else {
                  I0 = illum::BoundaryConditions(neigh_id);
                }
#else
                Type I0 = illum::BoundaryConditions(neigh_id);
#endif
                (*inter_coef)[face_block_id + num_out_face] = Vector3(I0, I0, I0); //значение на грани ( или коэффициенты интерполяции)
              }
              continue;
            }

            Vector3 I;
            // структура аналогичная  ::trace::GetLocNodes(...)
            for (ShortId num_node = 0; num_node < 3; ++num_node) {

              Vector3 &x = vec_x[num_cell].x[num_out_face][num_node];
              ShortId num_in_face = X0_ptr->in_face_id;

              Type I_x0 = GetIllumeFromInFace(neighbours[face_block_id + num_in_face], (*inter_coef)[face_block_id + num_in_face]);
              I[num_node] = GetIllum(x, X0_ptr->s, I_x0, grid.scattering[num_direction * count_cells + num_cell], *cell);

              X0_ptr++;
            } // num_node

            //записываем значение коэффициентов на грани
            (*inter_coef)[face_block_id + num_out_face] = I; // coef

            // если эта не граница, переносим это же значение входную грань соседней ячейки
            if (neigh_id >= 0) {
              (*inter_coef)[neigh_id] = I;
            }
          } // num_out_face
        }
        /*---------------------------------- конец FOR по ячейкам----------------------------------*/

        loc_norm = ReCalcIllum(num_direction, *inter_coef, grid, disp_illum[myid]);
      }
      /*--------------------------------конец первой секции---------------------------------*/

#pragma omp flush

#pragma omp single // nowait
      {
        MPI_Startall(section_1.requests_send.size(), section_1.requests_send.data());
        cuda::interface::CudaSendIllumAsync(section_1.size * n_illum, disp_illum[myid], grid.Illum);

        section_2.flags_send_to_gpu.assign(section_2.flags_send_to_gpu.size(), 0); // nowait
        MPI_Startall(section_2.requests_rcv.size(), section_2.requests_rcv.data());
        cuda::interface::CudaSyncStream(cuda::e_cuda_scattering_2);
      }

#pragma omp for
      for (IdType num_direction = section_1.size; num_direction < local_size; num_direction++) {

        const cell_local *X0_ptr = vec_x0[num_direction].data(); ///< индексация по массиву определяющих гранях (конвеерная т.к. заранее не известны позиции точек)
        const IntId *code_bound = inner_bound_code[num_direction].data();

        /*---------------------------------- далее FOR по ячейкам----------------------------------*/
        for (IdType h = 0; h < count_cells; ++h) {

          const IdType num_cell = sorted_id_cell[num_direction][h];
          const IdType face_block_id = num_cell * CELL_SIZE;

          elem_t *cell = &grid.cells[num_cell];

          // расчитываем излучения на выходящих гранях
          for (ShortId num_out_face = 0; num_out_face < CELL_SIZE; ++num_out_face) {

            const IdType neigh_id = neighbours[face_block_id + num_out_face]; ///< сосед к текущей грани

            // если эта грань входящая и граничная, то пропускаем её
            if (CHECK_BIT(face_states[num_direction][num_cell], num_out_face) == e_face_type_in) {
              if (neigh_id < 0) {
#ifdef USE_TRACE_THROUGH_INNER_BOUNDARY
                Type I0;
                if (neigh_id == e_bound_inner_source) {
                  Vector3 I_def = *code_bound >= 0 ? (*inter_coef)[*code_bound] : Vector3::Zero(); //т.е. определять будет ячейка, а не геометрия
                  I0 = illum::BoundaryConditions(neigh_id, *code_bound, I_def);
                  code_bound++;
                } else {
                  I0 = illum::BoundaryConditions(neigh_id);
                }
#else
                Type I0 = illum::BoundaryConditions(neigh_id);
#endif
                (*inter_coef)[face_block_id + num_out_face] = Vector3(I0, I0, I0); //значение на грани ( или коэффициенты интерполяции)
              }
              continue;
            }

            Vector3 I;
            // структура аналогичная  ::trace::GetLocNodes(...)
            for (ShortId num_node = 0; num_node < 3; ++num_node) {

              Vector3 &x = vec_x[num_cell].x[num_out_face][num_node];
              ShortId num_in_face = X0_ptr->in_face_id;

              Type I_x0 = GetIllumeFromInFace(neighbours[face_block_id + num_in_face], (*inter_coef)[face_block_id + num_in_face]);
              I[num_node] = GetIllum(x, X0_ptr->s, I_x0, grid.scattering[num_direction * count_cells + num_cell], *cell);

              X0_ptr++;
            } // num_node

            //записываем значение коэффициентов на грани
            (*inter_coef)[face_block_id + num_out_face] = I; // coef

            // если эта не граница, переносим это же значение входную грань соседней ячейки
            if (neigh_id >= 0) {
              (*inter_coef)[neigh_id] = I;
            }
          } // num_out_face
        }
        /*---------------------------------- конец FOR по ячейкам----------------------------------*/

        loc_norm = ReCalcIllum(num_direction, *inter_coef, grid, disp_illum[myid]);

#pragma omp critical
        {
          if (num_th == 0) //вместо критической секции пусть всегда первая нить управляет отправкой
          {
            // пересылаем первую пачку сообщений
            for (int i = 0; i < section_1.requests_rcv.size(); i++) {
              if (!section_1.flags_send_to_gpu[i]) {
                MPI_Test(&section_1.requests_rcv[i], &section_1.flags_send_to_gpu[i], &section_1.status_rcv[i]); //проверяем все запросы принятия сообщения
                if (section_1.flags_send_to_gpu[i])                                                              // если обмен завершён, но отправки не было
                {
                  const int src = section_1.status_rcv[i].MPI_TAG;
                  cuda::interface::CudaSendIllumAsync(n_illum * section_1.size, disp_illum[src], grid.Illum);
                }
              }
            }
          }

          MPI_Startall(np - 1, section_2.requests_send.data() + ((num_direction - section_1.size) * (np - 1)));

          if ((num_th == (count_th - 1)) && np > 1) // отличный от нулевого поток
          {
            for (int i = 0; i < section_2.requests_rcv.size(); i++) {
              if (!section_2.flags_send_to_gpu[i]) {
                MPI_Test(&section_2.requests_rcv[i], &section_2.flags_send_to_gpu[i], &section_2.status_rcv[i]); //проверяем все запросы принятия сообщения
                if (section_2.flags_send_to_gpu[i])                                                              // если обмен завершён, но отправки не было
                {
                  cuda::interface::CudaSendIllumAsync(n_illum, n_illum * (section_2.status_rcv[i].MPI_TAG), grid.Illum); //переслать данные на gpu
                }
              }
            }
          }
          cuda::interface::CudaSendIllumAsync(n_illum, ((local_disp + num_direction) * n_illum), grid.Illum);
        }
      }
      /*---------------------------------- конец FOR по направлениям----------------------------------*/

      if (loc_norm > norm) {
#pragma omp critical
        {
          if (loc_norm > norm) {
            norm = loc_norm;
          }
        }
      }
    } // parallel

    MPI_Request rq_norm;
    MPI_Iallgather(&norm, 1, MPI_DOUBLE, norms.data(), 1, MPI_DOUBLE, MPI_COMM_ILLUM, &rq_norm);

    bool ready = true;
    do {
      ready = true;
      for (int i = 0; i < section_1.requests_rcv.size(); i++) {
        if (!section_1.flags_send_to_gpu[i]) {
          ready = false;
          MPI_Test(&section_1.requests_rcv[i], &section_1.flags_send_to_gpu[i], &section_1.status_rcv[i]); //проверяем все запросы принятия сообщения
          if (section_1.flags_send_to_gpu[i])                                                              // если обмен завершён, но отправки не было
          {
            int src = (section_1.status_rcv[i].MPI_TAG);
            DIE_IF(src >= disp_illum.size())

            cuda::interface::CudaSendIllumAsync(n_illum * section_1.size, disp_illum[src], grid.Illum);
          }
        }
      }
    } while (!ready);

    do {
      ready = true;
      for (int i = 0; i < section_2.requests_rcv.size(); i++) {
        if (!section_2.flags_send_to_gpu[i]) {
          ready = false;
          MPI_Test(&section_2.requests_rcv[i], &section_2.flags_send_to_gpu[i], &section_2.status_rcv[i]); //проверяем все запросы принятия сообщения
          if (section_2.flags_send_to_gpu[i])                                                              // если обмен завершён, но отправки не было
          {
            cuda::interface::CudaSendIllumAsync(n_illum, n_illum * (section_2.status_rcv[i].MPI_TAG), grid.Illum); //переслать данные на gpu
          }
        }
      }

    } while (!ready);

    if (_solve_mode.max_number_of_iter >= 1) // пропуск первой итерации
    {
      cuda::interface::CalculateIntScatteringAsync(grid_direction, grid, 0, section_1.size, cuda::e_cuda_scattering_1);
      cuda::interface::CalculateIntScatteringAsync(grid_direction, grid, section_1.size, local_size, cuda::e_cuda_scattering_2);
    }

    MPI_Wait(&rq_norm, MPI_STATUS_IGNORE);
    for (auto n : norms)
      norm = std::max(norm, n);

    WRITE_LOG("End iter number#: %d, norm=%.16lf, time= %lf\n", iter, norm,
              (double)tick::duration_cast<tick::milliseconds>(tick::steady_clock::now() - start_clock).count() / 1000.);
    iter++;
  } while (norm > _solve_mode.accuracy && iter < _solve_mode.max_number_of_iter);

  return e_completion_success;
}
#endif //! defined ILLUM && defined SOLVERS  && !defined USE_MPI