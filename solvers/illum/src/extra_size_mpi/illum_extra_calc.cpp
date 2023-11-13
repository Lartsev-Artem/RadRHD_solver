/// \todo ReCalcIllum должен пересылать в локальный массив тот же что и ReCalcIllumLocal т.к.

#if defined ILLUM && defined SOLVERS && defined USE_MPI && defined USE_CUDA
#include "illum_extra_main.h"
#include "illum_mpi_sender.h"

#include "illum_utils.h"

#include "cuda_interface.h"

#include <omp.h>

#include <chrono>
namespace tick = std::chrono;

int illum::extra_size::CalculateIllum(const grid_directions_t &grid_direction, const std::vector<std::vector<bits_flag_t>> &face_states,
                                      const std::vector<IntId> &neighbours, const std::vector<std::vector<IntId>> &inner_bound_code,
                                      const std::vector<std::vector<cell_local>> &vec_x0, std::vector<BasePointTetra> &vec_x,
                                      const std::vector<std::vector<IntId>> &sorted_id_cell, grid_t &grid) {

  const IdType n_illum = grid.size;

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
    cuda::interface::CudaSyncStream(cuda::e_cuda_scattering_1);

    // if (myid == 0) {
    //   //самый высоких приоритет, т.к. надо расчитать, до конфликта с асинхронной отправкой
    //   cuda::interface::CalculateAllParamAsync(grid_direction, grid, cuda::e_cuda_params); //запустим расчёт параметров здесь
    //                                                                                       // на выходе получим ответ за 1 шаг до сходимости, но зато без ожидания на выходе
    // }

    /*---------------------------------- далее FOR по направлениям----------------------------------*/
    const IdType count_directions = grid_direction.size;

#pragma omp parallel default(none) firstprivate(count_directions, n_illum, myid, np, local_disp, local_size) \
    shared(sorted_id_cell, neighbours, face_states, vec_x0, vec_x, grid, norm, disp_illum, section_1,        \
           inner_bound_code)
    {

      const int count_th = omp_get_num_threads();
      const int num_th = omp_get_thread_num();

      const int count_cells = grid.size;

      Type loc_norm = -1;
      std::vector<Vector3> *inter_coef = &grid.inter_coef_all[omp_get_thread_num()]; ///< указатель на коэффициенты интерполяции по локальному для потока направлению

#pragma omp single // nowait
      {
        section_1.flags_send_to_gpu.assign(section_1.flags_send_to_gpu.size(), 0); // nowait
        MPI_Startall(section_1.requests_rcv.size(), section_1.requests_rcv.data());
        cuda::interface::CudaSyncStream(cuda::e_cuda_scattering_2);
      }

#pragma omp for
      for (IdType num_direction = 0; num_direction < local_size; num_direction++) {

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

        loc_norm = extra_size::ReCalcIllum(num_direction, *inter_coef, grid);

#pragma omp critical
        {
          MPI_Startall(np - 1, section_1.requests_send.data() + ((num_direction - 0) * (np - 1)));
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
    // MPI_Iallreduce(&norm, &norms[0], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_ILLUM, rq_norm);

    extra_size::ReCalcLocalIllum(grid_direction, grid);

    MPI_Waitall(section_1.requests_rcv.size(), section_1.requests_rcv.data(), MPI_STATUSES_IGNORE);

    if (_solve_mode.max_number_of_iter >= 1) // пропуск первой итерации
    {
      cuda::interface::CalculateIntScatteringMultiDev(grid_direction, grid, 0, 0, cuda::e_cuda_scattering_1);
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