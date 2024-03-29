#if defined SPECTRUM
#include "spectrum_full.h"
#ifdef SAVE_FULL_SPECTRUM

#include "illum_mpi_sender.h"

#include "cuda_interface.h"
#include "cuda_multi_interface.h"

#include <omp.h>

#include "global_consts.h"
#include "plunk.h"

#include "illum_utils.h"

namespace cuda_sep = cuda::interface::separate_device;

Type illum::full_spectrum::BoundaryConditions(const IdType type_bound, const Type frq0, const Type frq) {
  switch (type_bound) {

  case e_bound_free:
  case e_bound_lock:
    return 0;
  case e_bound_out_source:
#if GEOMETRY_TYPE == Cone
    return exp(B_Plank_log(1e12, frq, frq0) - LOG(kRadiation));
#endif
    return 0;

  case e_bound_inner_source:
#if GEOMETRY_TYPE == Sphere
    return exp(B_Plank_log(1e6, frq, frq0) - LOG(kRadiation));
#endif
    return 0;

  default:
    D_LD;
  }
}

Type illum::full_spectrum::ReCalcIllum(const IdType num_dir, const std::vector<std::vector<Type>> &inter_coef, grid_t &grid, const IdType dir_disp) {
  Type norm = -1;
  const IdType shift_dir = (grid.size_frq * (num_dir + dir_disp));
  const IdType shift_cell = grid.size_frq * grid.size_dir;
  IdType frq_idx = -1;
  IdType cell_idx = -1;

  for (IdType cell = 0; cell < grid.size; cell++) {
    IdType frq_id = cell * shift_cell + shift_dir;
    for (IdType frq = 0; frq < grid.size_frq; frq++) {

      Type curI = 0;
      for (size_t j = 0; j < CELL_SIZE; j++) {
        curI += inter_coef[grid.cells[cell].geo.id_faces[j]][frq]; //тут печаль с кэшами
      }
      curI /= CELL_SIZE;

      IdType id = frq_id + frq;
      norm = std::max(norm, fabs((grid.Illum[id] - curI) / curI));
      grid.Illum[id] = curI; //здесь по ячейкам
    }
  }

  // WRITE_LOG("max dF[%ld]: %lld\n", cell_idx, frq_idx);
  return norm;
}

int illum::full_spectrum::CalculateIllum(const grid_directions_t &grid_direction,
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
    const IdType count_frequencies = grid.size_frq;

#pragma omp parallel default(none) firstprivate(count_frequencies, count_directions, np, local_disp, local_size) \
    shared(sorted_graph, sorted_id_bound_face, inner_bound_code, vec_x0, grid, norm, section_1, grid_direction)
    {

      const int count_th = omp_get_num_threads();
      const int num_th = omp_get_thread_num();

      const IdType count_cells = grid.size;

      Type loc_norm = -1;
      std::vector<std::vector<Type>> *inter_coef = &grid.inter_coef_all[omp_get_thread_num()]; ///< указатель на коэффициенты интерполяции по локальному для потока направлению
#pragma omp single                                                                             // nowait
      {
        section_1.flags_send_to_gpu.assign(section_1.flags_send_to_gpu.size(), 0); // nowait
        MPI_Startall(section_1.requests_rcv.size(), section_1.requests_rcv.data());
        cuda::interface::CudaSyncStream(cuda::e_cuda_scattering_1);
      }

#pragma omp for
      for (IdType num_direction = 0; num_direction < local_size; num_direction++) {
        /*---------------------------------- FOR по граничным граням----------------------------------*/
        for (auto num_face : sorted_id_bound_face[num_direction]) {
          const IdType neigh_id = grid.faces[num_face].geo.id_r; //т.к. мы проходим границу здесь будет граничный признак

          const Type *frq = grid.frq_grid.data();
          std::vector<Type> *init_frq_coef = &(*inter_coef)[num_face];

          for (IdType num_frq = 0; num_frq < count_frequencies; num_frq++) {
            Type frq0 = *(frq);
            frq++;
            (*init_frq_coef)[num_frq] = full_spectrum::BoundaryConditions(neigh_id, frq0, *frq); //значение на грани ( или коэффициенты интерполяции)

            // if ((*init_frq_coef)[num_frq] != 0.0)
            // {
            //   WRITE_LOG("frq[%lf-%lf]= %lf\n", frq0, *frq, (*init_frq_coef)[num_frq]);
            // }
          }
        }
        // индексация по массиву определяющих гранях (конвейерная т.к. заранее не известны позиции точек)
        const Type *s = vec_x0[num_direction].s.data();
        const face_loc_id_t *in_face = vec_x0[num_direction].in_face_id.data();
        /*---------------------------------- далее FOR по неосвященным граням----------------------------------*/
        for (auto fc_pair : sorted_graph[num_direction]) {

          const IntId num_cell = fc_pair.cell;
          const uint32_t num_loc_face = fc_pair.loc_face;
          elem_t *cell = &grid.cells[num_cell];
          const Vector3 &x = cell->geo.center;
          const Type *S = &grid.scattering[count_frequencies * (num_cell * local_size + num_direction)];
          const face_loc_id_t id_in_faces = *(in_face);
          ++in_face;

          /*---------------------------------- далее FOR по частотам----------------------------------*/
          cell->cell_data->InitDirection(grid_direction.directions[num_direction].dir);
          std::vector<Type> *frq_coef_A = &(*inter_coef)[cell->geo.id_faces[id_in_faces.a]];
          std::vector<Type> *frq_coef_B = &(*inter_coef)[cell->geo.id_faces[id_in_faces.b]];
          std::vector<Type> *frq_coef_C = &(*inter_coef)[cell->geo.id_faces[id_in_faces.c]];
          std::vector<Type> *frq_coef_Res = &(*inter_coef)[cell->geo.id_faces[num_loc_face]];
          const Type *frq = grid.frq_grid.data();

          for (IdType num_frq = 0; num_frq < count_frequencies; num_frq++) {
            Type frq0 = *(frq);
            frq++;
            Type k;
            const Type rhs = GetRhsOpt(x, S[num_frq], *cell, k, frq0, *frq);

            alignas(32) Type I0[4] = {(*frq_coef_A)[num_frq], (*frq_coef_B)[num_frq], (*frq_coef_C)[num_frq], 0};
            (*frq_coef_Res)[num_frq] = GetIllum(I0, s, k, rhs);
          }
          /*---------------------------------- конец FOR по частотам----------------------------------*/
          s += NODE_SIZE;
        }
        /*---------------------------------- конец FOR по ячейкам----------------------------------*/
        loc_norm = ReCalcIllum(num_direction, *inter_coef, grid, local_disp);

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

    if (_solve_mode.max_number_of_iter >= 1) // пропуск первой итерации
    {
      cuda_sep::CalculateFullSpectrumIntScattering(grid_direction, grid, 0, local_size, cuda::e_cuda_scattering_1);
    }

    MPI_Wait(&rq_norm, MPI_STATUS_IGNORE);
    norm = fabs(glob_norm);

    WRITE_LOG("End iter number#: %d, norm=%.16lf, time= %lf\n", iter, norm, time.get_delta_time_sec());
    iter++;
  } while (norm > _solve_mode.accuracy && iter < _solve_mode.max_number_of_iter);

  return e_completion_success;
}
#endif //! SAVE_FULL_SPECTRUM
#endif //! SPECTRUM