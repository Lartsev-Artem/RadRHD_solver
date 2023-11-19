#if defined ILLUM && defined SOLVERS //&& !defined USE_MPI
#include "illum_calc_cpu.h"

#ifndef TRANSFER_CELL_TO_FACE
#include "illum_params.h"
#include "illum_utils.h"
#include "scattering.h"

#include "solvers_struct.h"

#ifdef USE_CUDA
#include "cuda_interface.h"
#endif

#include <omp.h>

#include <chrono>
namespace tick = std::chrono;

#ifdef INTERPOLATION_ON_FACES
static const BaseTetra_t tetra;
#endif

int illum::cpu::CalculateIllum(const grid_directions_t &grid_direction, const std::vector<std::vector<bits_flag_t>> &face_states,
                               const std::vector<IntId> &neighbours, const std::vector<std::vector<IntId>> &inner_bound_code,
                               const std::vector<std::vector<cell_local>> &vec_x0, std::vector<BasePointTetra> &vec_x,
                               const std::vector<std::vector<IntId>> &sorted_id_cell, grid_t &grid) {

  // вроде не обязательно. ПРОВЕРИТЬ
  // Illum.assign(4 * count_cells * count_directions, 0);
  // int_scattering.assign(count_cells * count_directions, 0);

  int iter = 0;    ///< номер итерации
  double norm = 0; ///< норма ошибки

  do {
    auto start_clock = tick::steady_clock::now();

    norm = -1;
    /*---------------------------------- далее FOR по направлениям----------------------------------*/
    const IdType count_directions = grid_direction.size;

#pragma omp parallel default(none) firstprivate(count_directions) shared(sorted_id_cell, neighbours, face_states, inner_bound_code, vec_x0, vec_x, grid, norm)
    {
      const IdType count_cells = grid.size;

      Type loc_norm = -1;
      std::vector<Vector3> *inter_coef = &grid.inter_coef_all[omp_get_thread_num()]; ///< указатель на коэффициенты интерполяции по локальному для потока направлению

#pragma omp for
      for (IdType num_direction = 0; num_direction < count_directions; ++num_direction) {

        const cell_local *X0_ptr = vec_x0[num_direction].data(); ///< индексация по массиву определяющих гранях (конвеерная т.к. заранее не известны позиции точек)
        const IntId *code_bound = inner_bound_code[num_direction].data();
        /*---------------------------------- далее FOR по ячейкам----------------------------------*/
        for (IdType h = 0; h < count_cells; ++h) {

          const IdType num_cell = sorted_id_cell[num_direction][h];
          elem_t *cell = &grid.cells[num_cell];

          // расчитываем излучения на выходящих гранях
          for (ShortId num_out_face = 0; num_out_face < CELL_SIZE; ++num_out_face) {

            const IdType neigh_id = neighbours[num_cell * CELL_SIZE + num_out_face]; ///< сосед к текущей грани

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
                (*inter_coef)[num_cell * CELL_SIZE + num_out_face] = Vector3(I0, I0, I0); //значение на грани ( или коэффициенты интерполяции)
              }
              continue;
            }

            Vector3 I;
            // структура аналогичная  ::trace::GetLocNodes(...)
            for (ShortId num_node = 0; num_node < 3; ++num_node) {

              Vector3 &x = vec_x[num_cell].x[num_out_face][num_node];
              ShortId num_in_face = X0_ptr->in_face_id;
#ifdef INTERPOLATION_ON_FACES
              Type I_x0 = GetIllumeFromInFace(neighbours[num_cell * CELL_SIZE + num_in_face], (*inter_coef)[num_cell * CELL_SIZE + num_in_face], X0_ptr->x0);
#else
              Type I_x0 = GetIllumeFromInFace(neighbours[num_cell * CELL_SIZE + num_in_face], (*inter_coef)[num_cell * CELL_SIZE + num_in_face]);
#endif
              I[num_node] = GetIllum(x, X0_ptr->s, I_x0, grid.scattering[num_direction * count_cells + num_cell], *cell);

              X0_ptr++;
            } // num_node

#ifdef INTERPOLATION_ON_FACES
            Type i3 = (I[0] + I[1] + I[2]) / 3;
            grid.Illum[CELL_SIZE * (num_direction * grid.size + num_cell) + num_out_face] = i3;
            if (neigh_id >= 0) {
              grid.Illum[CELL_SIZE * (num_direction * grid.size) + neigh_id] = i3;
            }
            // если хранить не значения у коэфф. интерполяции
            if (num_out_face == 3)
              I = tetra.inclined_face_inverse * I; // GetInterpolationCoefInverse(inclined_face_inverse, I);
            else
              I = tetra.straight_face_inverse * I; // GetInterpolationCoefInverse(straight_face_inverse, I);
#endif

            //записываем значение коэффициентов на грани
            (*inter_coef)[num_cell * CELL_SIZE + num_out_face] = I; // coef

            // если эта не граница, переносим это же значение входную грань соседней ячейки
            if (neigh_id >= 0) {
              (*inter_coef)[neigh_id] = I;
            }
          } // num_out_face
        }
        /*---------------------------------- конец FOR по ячейкам----------------------------------*/
#ifndef INTERPOLATION_ON_FACES
        loc_norm = ReCalcIllum(num_direction, *inter_coef, grid);
#endif
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

    if (_solve_mode.max_number_of_iter >= 1) // пропуск первой итерации
    {
#ifndef USE_CUDA
      scattering::CalculateIntCPU(grid_direction, grid);
#else
      cuda::interface::CalculateIntScattering(grid_direction, grid);
#endif
    }

    WRITE_LOG("End iter number#: %d, norm=%.16lf, time= %lf\n", iter, norm,
              (double)tick::duration_cast<tick::milliseconds>(tick::steady_clock::now() - start_clock).count() / 1000.);
    iter++;
  } while (norm > _solve_mode.accuracy && iter < _solve_mode.max_number_of_iter);

  return e_completion_success;
}

void illum::cpu::CalculateIllumParam(const grid_directions_t &grid_direction, grid_t &grid) {

#if !defined USE_CUDA
  GetEnergy(grid_direction, grid);
  GetStream(grid_direction, grid);
  GetImpuls(grid_direction, grid);
#else
  cuda::interface::CalculateAllParam(grid_direction, grid);
#endif
}

int illum::cpu::CalculateAdditionalIllum(const grid_directions_t &grid_direction,
                                         const std::vector<std::vector<bits_flag_t>> &face_states,
                                         const std::vector<IntId> &neighbours, const std::vector<std::vector<IntId>> &inner_bound_code,
                                         const std::vector<std::vector<cell_local>> &vec_x0, std::vector<BasePointTetra> &vec_x,
                                         const std::vector<std::vector<IntId>> &sorted_id_cell, grid_t &grid) {

  auto start_clock = tick::steady_clock::now();

  /*---------------------------------- далее FOR по направлениям----------------------------------*/

  const IdType size_directions = grid_direction.loc_size;

#pragma omp parallel default(none) firstprivate(size_directions) shared(sorted_id_cell, neighbours, face_states, inner_bound_code, vec_x0, vec_x, grid)
  {
    const IdType count_cells = grid.size;

    Type loc_norm = -1;
    std::vector<Vector3> *inter_coef = &grid.inter_coef_all[omp_get_thread_num()]; ///< указатель на коэффициенты интерполяции по локальному для потока направлению

#pragma omp for
    for (IdType num_direction = 0; num_direction < size_directions; ++num_direction) {

      const cell_local *X0_ptr = vec_x0[num_direction].data(); ///< индексация по массиву определяющих гранях (конвеерная т.к. заранее не известны позиции точек)
      const IntId *code_bound = inner_bound_code[num_direction].data();
      /*---------------------------------- далее FOR по ячейкам----------------------------------*/
      for (IdType h = 0; h < count_cells; ++h) {

        const IdType num_cell = sorted_id_cell[num_direction][h];
        elem_t *cell = &grid.cells[num_cell];

        // расчитываем излучения на выходящих гранях
        for (ShortId num_out_face = 0; num_out_face < CELL_SIZE; ++num_out_face) {

          const IdType neigh_id = neighbours[num_cell * CELL_SIZE + num_out_face]; ///< сосед к текущей грани

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
              (*inter_coef)[num_cell * CELL_SIZE + num_out_face] = Vector3(I0, I0, I0); //значение на грани ( или коэффициенты интерполяции)
            }
            continue;
          }

          Vector3 I;
          // структура аналогичная  ::trace::GetLocNodes(...)
          for (ShortId num_node = 0; num_node < 3; ++num_node) {

            Vector3 &x = vec_x[num_cell].x[num_out_face][num_node];
            ShortId num_in_face = X0_ptr->in_face_id;

            Type I_x0 = GetIllumeFromInFace(neighbours[num_cell * CELL_SIZE + num_in_face], (*inter_coef)[num_cell * CELL_SIZE + num_in_face]);

            I[num_node] = GetIllum(x, X0_ptr->s, I_x0, grid.scattering[num_direction * count_cells + num_cell], *cell);

            X0_ptr++;
          } // num_node

          //записываем значение коэффициентов на грани
          (*inter_coef)[num_cell * CELL_SIZE + num_out_face] = I; // coef

          // если эта не граница, переносим это же значение входную грань соседней ячейки
          if (neigh_id >= 0) {
            (*inter_coef)[neigh_id] = I;
          }
        } // num_out_face
      }
      /*---------------------------------- конец FOR по ячейкам----------------------------------*/

      loc_norm = ReCalcIllum(num_direction, *inter_coef, grid);
    }
    /*---------------------------------- конец FOR по направлениям----------------------------------*/

  } // parallel

  WRITE_LOG("End iter number#:time= %lf\n",
            (double)tick::duration_cast<tick::milliseconds>(tick::steady_clock::now() - start_clock).count() / 1000.);

  return e_completion_success;
}

#endif
#endif //! defined ILLUM && defined SOLVERS  && !defined USE_MPI