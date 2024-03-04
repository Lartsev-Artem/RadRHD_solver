#if defined ILLUM && defined SOLVERS && defined USE_MPI && defined USE_CUDA

#include "illum_calc_gpu_async.h"
#if defined TRANSFER_CELL_TO_FACE && defined SEPARATE_GPU
#include "global_consts.h"
#include "illum_mpi_sender.h"
#include "illum_utils.h"
#include "plunk.h"

#include "cuda_interface.h"
#include "cuda_multi_interface.h"

#include <omp.h>

#include <chrono>
namespace tick = std::chrono;
namespace cuda_sep = cuda::interface::separate_device;

#define builtin_prefetch(...) __builtin_prefetch(__VA_ARGS__)
#include "spectrum_utils.h"
int illum::separate_gpu::CalculateIllum(const grid_directions_t &grid_direction,
                                        const std::vector<std::vector<IntId>> &inner_bound_code,
                                        const std::vector<align_cell_local> &vec_x0,
                                        const std::vector<std::vector<graph_pair_t>> &sorted_graph,
                                        const std::vector<std::vector<IntId>> &sorted_id_bound_face,
                                        grid_t &grid) {

  int np = get_mpi_np();
  int myid = get_mpi_id();

  DIE_IF(np <= 1); //на одном узле не работает. Тогда надо закрывать пересылки mpi

  const IdType local_size = grid_direction.loc_size;
  const IdType local_disp = grid_direction.loc_shift;

#ifdef SPECTRUM
  const Type frq0 = grid.frq_grid[grid.cur_frq_id];
  const Type frq1 = grid.frq_grid[grid.cur_frq_id + 1];
#endif

  int iter = 0;    ///< номер итерации
  double norm = 0; ///< норма ошибки

  do {
    Timer time;
    time.start_timer();
    norm = -1;

    if (section_1.requests_send[0] != MPI_REQUEST_NULL) {
      MPI_Waitall(section_1.requests_send.size(), section_1.requests_send.data(), MPI_STATUSES_IGNORE);
    }
    cuda::interface::CudaSyncStream(cuda::e_cuda_scattering_1);

#ifdef MULTI_GPU //если разделение в порядке очереди надо запускать расчёт до пересылки новых данных
    if (myid == 0) {
      //самый высоких приоритет, т.к. надо расчитать, до конфликта с асинхронной отправкой
      cuda::interface::CalculateAllParamAsync(grid_direction, grid, cuda::e_cuda_params); //запустим расчёт параметров здесь
                                                                                          // на выходе получим ответ за 1 шаг до сходимости, но зато без ожидания на выходе
    }
#endif

    /*---------------------------------- далее FOR по направлениям----------------------------------*/
    const IdType count_directions = grid_direction.size;

#ifdef SPECTRUM
#pragma omp parallel default(none) firstprivate(count_directions, np, local_disp, local_size, frq0, frq1) \
    shared(sorted_graph, sorted_id_bound_face, inner_bound_code, vec_x0, grid, norm, section_1, grid_direction, glb_files, log_enable)
#else
#pragma omp parallel default(none) firstprivate(count_directions, np, local_disp, local_size) \
    shared(sorted_graph, sorted_id_bound_face, inner_bound_code, vec_x0, grid, norm, section_1)
#endif
    {
      const int count_th = omp_get_num_threads();
      const int num_th = omp_get_thread_num();

      const IdType count_cells = grid.size;

      Type loc_norm = -1;
      std::vector<Type> *inter_coef = &grid.inter_coef_all[omp_get_thread_num()]; ///< указатель на коэффициенты интерполяции по локальному для потока направлению
      alignas(32) Type I0[4];

#pragma omp single // nowait
      {
        section_1.flags_send_to_gpu.assign(section_1.flags_send_to_gpu.size(), 0); // nowait
        MPI_Startall(section_1.requests_rcv.size(), section_1.requests_rcv.data());
        cuda::interface::CudaSyncStream(cuda::e_cuda_scattering_1);
      }

#pragma omp for
      for (IdType num_direction = 0; num_direction < local_size; num_direction++) {
        /*---------------------------------- FOR по граничным граням----------------------------------*/
#ifdef USE_TRACE_THROUGH_INNER_BOUNDARY
        const IntId *code_bound = inner_bound_code[num_direction].data();
#pragma error "need recalc code bound from cell to face"
#endif
        for (auto num_face : sorted_id_bound_face[num_direction]) {
          const IdType neigh_id = grid.faces[num_face].geo.id_r;         //т.к. мы проходим границу здесь будет граничный признак
          (*inter_coef)[num_face] = illum::BoundaryConditions(neigh_id); //значение на грани ( или коэффициенты интерполяции)
        }
#ifdef LOG_SPECTRUM
        if (num_direction == 32 + local_disp)
          log_enable = 1;
#endif

        // индексация по массиву определяющих гранях (конвеерная т.к. заранее не известны позиции точек)
        const Type *s = vec_x0[num_direction].s.data();
        const face_loc_id_t *in_face = vec_x0[num_direction].in_face_id.data();
        /*---------------------------------- далее FOR по неосвященным граням----------------------------------*/
        // for (auto fc_pair : sorted_graph[num_direction])
        const graph_pair_t *fc_pair = sorted_graph[num_direction].data();
        const IdType Ncells = sorted_graph[num_direction].size();
        for (IdType i = 0; i < Ncells; i++) {

          const IntId num_cell = fc_pair->cell;
          const uint32_t num_loc_face = fc_pair->loc_face;

          builtin_prefetch(&(grid.scattering[num_cell * local_size + num_direction]), 0, 0);

          elem_t *cell = &grid.cells[num_cell];
          const Vector3 &x = cell->geo.center; // vec_x[num_cell].x[num_loc_face][num_node];
#ifdef SPECTRUM
          cell->cell_data->InitDirection(grid_direction.directions[num_direction].dir);
#endif
          const face_loc_id_t id_in_faces = *(in_face);
          ++in_face;

          builtin_prefetch(&((*inter_coef)[cell->geo.id_faces[id_in_faces.a]]), 0, 0);
          builtin_prefetch(&((*inter_coef)[cell->geo.id_faces[id_in_faces.b]]), 0, 0);
          builtin_prefetch(&((*inter_coef)[cell->geo.id_faces[id_in_faces.c]]), 0, 0);
          builtin_prefetch(I0, 1, 1);
          Type k;
          const Type S = grid.scattering[num_cell * local_size + num_direction];
#ifdef SPECTRUM
          const Type rhs = GetRhsOpt(x, S, *cell, k, frq0, frq1);
#else
          const Type rhs = GetRhsOpt(x, S, *cell, k);
#endif

          I0[0] = (*inter_coef)[cell->geo.id_faces[id_in_faces.a]];
          I0[1] = (*inter_coef)[cell->geo.id_faces[id_in_faces.b]];
          I0[2] = (*inter_coef)[cell->geo.id_faces[id_in_faces.c]];

          builtin_prefetch(&((*inter_coef)[cell->geo.id_faces[num_loc_face]]), 1, 0);

          ++fc_pair;
          builtin_prefetch(&(grid.cells[fc_pair->cell]), 1, 2); //грузим следующую ячейку

          Type I;
          if (k > numeric_limit_abs_coef) {
            I = GetIllum(I0, s, k, rhs);
          } else {
            I = GetIllumLimit(I0, s, k, rhs);
          }
#ifdef DEBUG
          if (std::isnan(fabs(I)) || std::isinf(fabs(I)) || I < 0) {
            WRITE_LOG_ERR("nan[%d %d], k=%e, rhs=%e, I0:%e, %e, %e, %e\n", num_direction, num_cell, k, rhs, I0[0], I0[1], I0[2], S);
            D_LD;
          }
#endif
          (*inter_coef)[cell->geo.id_faces[num_loc_face]] = I;
          s += NODE_SIZE;
        }
        /*---------------------------------- конец FOR по ячейкам----------------------------------*/

        loc_norm = separate_gpu::ReCalcIllumOpt(num_direction, *inter_coef, grid, local_disp);

#pragma omp critical
        {
          MPI_Startall(np - 1, section_1.requests_send.data() + ((num_direction - 0) * (np - 1)));

          if (loc_norm > norm) {
            norm = loc_norm;
          }
        }
      }
      /*---------------------------------- конец FOR по направлениям----------------------------------*/
    } // parallel

    MPI_Request rq_norm;
    Type glob_norm = -1;
    MPI_Iallreduce(&norm, &glob_norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_ILLUM, &rq_norm);

    MPI_Waitall(section_1.requests_rcv.size(), section_1.requests_rcv.data(), MPI_STATUSES_IGNORE);

    if (LIKELY(_solve_mode.max_number_of_iter >= 1)) // пропуск первой итерации
    {
#ifdef SPECTRUM
      cuda_sep::CalculateSpectrumIntScattering(grid_direction, grid, 0, local_size, cuda::e_cuda_scattering_1);
#else
      cuda_sep::CalculateIntScatteringAsync(grid_direction, grid, 0, local_size, cuda::e_cuda_scattering_1);
#endif
    }

    MPI_Wait(&rq_norm, MPI_STATUS_IGNORE);
    norm = fabs(glob_norm);

    // WRITE_LOG("End iter number#: %d, norm=%.16lf, time= %lf\n", iter, norm, time.get_delta_time_sec());
    iter++;
  } while (norm > _solve_mode.accuracy && iter < _solve_mode.max_number_of_iter);

  return e_completion_success;
}
#endif //! TRANSFER_CELL_TO_FACE
#endif //! defined ILLUM && defined SOLVERS  && !defined USE_MPI