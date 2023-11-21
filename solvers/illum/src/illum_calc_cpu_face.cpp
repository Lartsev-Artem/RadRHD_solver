#if defined ILLUM && defined SOLVERS //&& !defined USE_MPI
#include "illum_calc_cpu.h"

#ifdef TRANSFER_CELL_TO_FACE
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
#include "reader_bin.h"
int illum::cpu::CalculateIllumFace(const grid_directions_t &grid_direction,
                                   const std::vector<std::vector<IntId>> &inner_bound_code,
                                   const std::vector<std::vector<cell_local>> &vec_x0, std::vector<BasePointTetra> &vec_x,
                                   const std::vector<std::vector<graph_pair_t>> &sorted_graph,
                                   const std::vector<std::vector<IntId>> &sorted_id_bound_face,
                                   grid_t &grid) {

  int iter = 0;    ///< номер итерации
  double norm = 0; ///< норма ошибки

  do {
    auto start_clock = tick::steady_clock::now();

    norm = -1;
    /*---------------------------------- далее FOR по направлениям----------------------------------*/
    const IdType count_directions = grid_direction.size;

#pragma omp parallel default(none) firstprivate(count_directions) shared(sorted_graph, sorted_id_bound_face, inner_bound_code, vec_x0, vec_x, grid, norm)
    {
      const IdType count_cells = grid.size;

      Type loc_norm = -1;
      std::vector<Type> *inter_coef = &grid.inter_coef_all[omp_get_thread_num()]; ///< указатель на коэффициенты интерполяции по локальному для потока направлению

#pragma omp for
      for (IdType num_direction = 0; num_direction < count_directions; ++num_direction) {

        (*inter_coef).assign(inter_coef->size(), -100);
        /*---------------------------------- FOR по граничным граням----------------------------------*/
#ifdef USE_TRACE_THROUGH_INNER_BOUNDARY
        const IntId *code_bound = inner_bound_code[num_direction].data();
#pragma error "need recalc code bound from cell to face"
#endif
        for (auto num_face : sorted_id_bound_face[num_direction]) {
          const IdType neigh_id = grid.faces[num_face].geo.id_r;         //т.к. мы проходим границу здесь будет граничный признак
          (*inter_coef)[num_face] = illum::BoundaryConditions(neigh_id); //значение на грани ( или коэффициенты интерполяции)
        }

        const cell_local *X0_ptr = vec_x0[num_direction].data(); ///< индексация по массиву определяющих гранях (конвеерная т.к. заранее не известны позиции точек)
        /*---------------------------------- далее FOR по неосвященным граням----------------------------------*/
        for (auto &fc_pair : sorted_graph[num_direction]) {

          IntId num_cell = fc_pair.cell;
          uint32_t num_loc_face = fc_pair.loc_face;
          elem_t *cell = &grid.cells[num_cell];

          Vector3 I;
          // структура аналогичная  ::trace::GetLocNodes(...)
          for (ShortId num_node = 0; num_node < 3; ++num_node) {

            IntId in_face_id = cell->geo.id_faces[X0_ptr->in_face_id]; //из локального номера в глобальный

#ifdef INTERPOLATION_ON_FACES
#pragma error "need convert coef interpolation to face_value"
#else
            Type I_x0 = (*inter_coef)[in_face_id]; //Здесь даже если попадаем на границу, она должна быть определена.
#endif
            const Vector3 &x = vec_x[num_cell].x[num_loc_face][num_node];
            I[num_node] = GetIllum(x, X0_ptr->s, I_x0, grid.scattering[num_direction * count_cells + num_cell], *cell);
            X0_ptr++;
          } // num_node

          //записываем значение коэффициентов на грани
          (*inter_coef)[cell->geo.id_faces[num_loc_face]] = (I[0] + I[1] + I[2]) / 3.;
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

#endif // calc on face
#endif //! defined ILLUM && defined SOLVERS  && !defined USE_MPI
