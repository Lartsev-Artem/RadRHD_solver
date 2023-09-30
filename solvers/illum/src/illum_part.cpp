#if defined ILLUM && defined SOLVERS // && !defined USE_MPI
#include "illum_part.h"
#include "illum_utils.h"
#include "scattering.h"

#include "geo_types.h"

#include <omp.h>

#include <chrono>
namespace tick = std::chrono;

static std::vector<std::vector<Vector3>> inter_coef_all; ///< коэффициенты интерполяции локальные для каждого потока

static Type CalculateIllumeOnInnerFace(const int num_in_face, const std::vector<face_t> &faces, elem_t *cell, Vector3 &inter_coef) {
  const int id_face_ = cell->geo.id_faces[num_in_face]; //номер грани
  const int neigh_id = faces[id_face_].geo.id_r;        // признак ГУ + связь с глобальной нумерацией

  if (neigh_id < 0) {
    Type I_x0 = illum::BoundaryConditions(neigh_id);
    inter_coef = Vector3(I_x0, I_x0, I_x0);
    return I_x0;
  }
  // Vector3 coef = grid[num_cell].nodes_value[num_in_face];
  Vector3 coef = inter_coef; // cell->illum_val.coef_inter[num_in_face];

  /// \note сейчас храним значения а не коэффициента интерполяции

  // Vector2	x0_local = X0[ShiftX0 + posX0++]; // grid[num_cell].x0_loc[num_in_face_dir];
  // I_x0 = x0_local[0] * coef[0] + x0_local[1] * coef[1] + coef[2];

  Type I_x0 = (coef[0] + coef[1] + coef[2]) / 3.;

  if (I_x0 < 0) {
    D_L;
    return 0;
  }
  return I_x0;
}

int illum::CalculateIllum(const grid_directions_t &grid_direction, const std::vector<std::vector<bits_flag_t>> &face_states,
                          const std::vector<IntId> &neighbours,
                          const std::vector<std::vector<cell_local>> &vec_x0, std::vector<BasePointTetra> &vec_x,
                          const std::vector<std::vector<IntId>> &sorted_id_cell, grid_t &grid) {

  // вроде не обязательно. ПРОВЕРИТЬ
  // Illum.assign(4 * count_cells * count_directions, 0);
  // int_scattering.assign(count_cells * count_directions, 0);

  int iter = 0;    ///< номер итерации
  double norm = 0; ///< норма ошибки

  inter_coef_all.resize(omp_get_max_threads());
  for (size_t i = 0; i < inter_coef_all.size(); i++) {
    inter_coef_all[i].resize(grid.size * CELL_SIZE);
  }

  do {
    auto start_clock = tick::steady_clock::now();

    norm = -1;
    /*---------------------------------- далее FOR по направлениям----------------------------------*/

#pragma omp parallel default(none) shared(sorted_id_cell, neighbours, face_states, vec_x0, vec_x, grid, norm, inter_coef_all)
    {
      const int count_directions = grid_direction.size;
      const int count_cells = grid.size;

      Type loc_norm = -1;
      std::vector<Vector3> *inter_coef = &inter_coef_all[omp_get_thread_num()]; ///< указатель на коэффициенты интерполяции по локальному для потока направлению

#pragma omp for
      for (int num_direction = 0; num_direction < count_directions; ++num_direction) {

        const cell_local *X0_ptr = vec_x0[num_direction].data(); ///< индексация по массиву определяющих гранях (конвеерная т.к. заранее не известны позиции точек)

        /*---------------------------------- далее FOR по ячейкам----------------------------------*/
        for (int h = 0; h < count_cells; ++h) {

          const int num_cell = sorted_id_cell[num_direction][h];
          elem_t *cell = &grid.cells[num_cell];

          // расчитываем излучения на выходящих гранях
          for (ShortId num_out_face = 0; num_out_face < CELL_SIZE; ++num_out_face) {

            // если эта грань входящая и не граничная, то пропускаем её
            if (CHECK_BIT(face_states[num_direction][num_cell], num_out_face) == e_face_type_in) {
              if (neighbours[num_cell * CELL_SIZE + num_out_face] < 0) {
                Type I0 = illum::BoundaryConditions(neighbours[num_cell * CELL_SIZE + num_out_face]);
                (*inter_coef)[num_cell * CELL_SIZE + num_out_face] = Vector3(I0, I0, I0); //значение на грани ( или коэффициенты интерполяции)
              }
              continue;
            }

            Vector3 I;
            // структура аналогичная  ::trace::GetLocNodes(...)
            for (int num_node = 0; num_node < 3; ++num_node) {

              Vector3 &x = vec_x[num_cell].x[num_out_face][num_node];
              ShortId num_in_face = X0_ptr->in_face_id;

              Type I_x0 = CalculateIllumeOnInnerFace(num_in_face, grid.faces, cell, (*inter_coef)[num_cell * CELL_SIZE + num_in_face]);
              I[num_node] = GetIllum(x, X0_ptr->s, I_x0, grid.scattering[num_direction * count_cells + num_cell], *cell);

              X0_ptr++;
            } // num_node

            // если хранить не значения у коэфф. интерполяции
            {
                // Vector3 coef;
                // if (num_out_face == 3)
                //	coef = inclined_face_inverse * I;// GetInterpolationCoefInverse(inclined_face_inverse, I);
                // else
                //	coef = straight_face_inverse * I;// GetInterpolationCoefInverse(straight_face_inverse, I);
            }

            //записываем значение коэффициентов на грани
            (*inter_coef)[num_cell * CELL_SIZE + num_out_face] = I; // coef

            // если эта не граница, переносим это же значение входную грань соседней ячейки
            const int id_face = neighbours[num_cell * CELL_SIZE + num_out_face];
            if (id_face >= 0) {
              (*inter_coef)[id_face] = I;
            }
          } // num_out_face
        }
        /*---------------------------------- конец FOR по ячейкам----------------------------------*/

        loc_norm = ReCalcIllum(num_direction, *inter_coef, grid);
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

    if (iter > 0) // пропуск первой итерации
    {
      scattering::CalculateIntCPU(grid_direction, grid);
    }

    WRITE_LOG("End iter number#: %d, norm=%lf, time= %lf\n", iter, norm,
              (double)tick::duration_cast<tick::milliseconds>(tick::steady_clock::now() - start_clock).count() / 1000.);
    iter++;
  } while (norm > _solve_mode.accuracy && iter < _solve_mode.max_number_of_iter);

  return e_completion_success;
}

#endif //! defined ILLUM && defined SOLVERS  && !defined USE_MPI