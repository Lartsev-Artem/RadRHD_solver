#include "solvers_config.h"
#ifdef SPECTRUM
#include "spec_all.h"

#ifdef USE_CUDA
#include "cuda_interface.h"
#endif
#include <omp.h>

#define USE_ONE_DIRECTION

#include "global_value.h"
#include "plunk.h"

#include "writer_txt.h"
static void WriteRes(const std::string &name_file, std::vector<Type> &y) {
  files_sys::txt::WriteSimple(name_file, y);
}

static void WriteRes(const std::string &name_file, std::vector<Type> &x, std::vector<Type> &y) {
  files_sys::txt::WriteSimple(name_file, x, y);
}

static Type get_full_illum(IdType count_dir, IdType num_frq, const grid_t &grid) {
  const IdType N = grid.size;
  const IdType shift_dir = 0;
  Type sumI = 0;
  for (size_t cell = 0; cell < N; cell++) {
    for (size_t j = 0; j < CELL_SIZE; j++) {
      sumI += grid.Illum[CELL_SIZE * (shift_dir + cell) + j];
    }
  }
  return sumI / kDist;
}

int illum::spec::CalculateIllumFace(const grid_directions_t &grid_direction,
                                    const std::vector<std::vector<IntId>> &inner_bound_code,
                                    const std::vector<align_cell_local> &vec_x0,
                                    const std::vector<std::vector<graph_pair_t>> &sorted_graph,
                                    const std::vector<std::vector<IntId>> &sorted_id_bound_face,
                                    grid_t &grid) {

  int iter = 0;    ///< номер итерации
  double norm = 0; ///< норма ошибки

  std::vector<Type> spectrum;
  get_splitting_spectrum(spectrum);

  const IdType count_frequencies = spectrum.size(); ///< число разбиений по частоте

  std::vector<Type> res_spectrum(count_frequencies - 1);

  for (IdType num_frq = 0; num_frq < count_frequencies - 1; num_frq++) {

    Type frq0 = spectrum[num_frq];
    Type frq1 = spectrum[num_frq + 1];
    // WRITE_LOG("frq= %lf %lf\n", frq0, frq1);
    do {

      norm = -1;
      /*---------------------------------- далее FOR по направлениям----------------------------------*/
#ifdef USE_ONE_DIRECTION
      const IdType count_directions = 1;
#else
      const IdType count_directions = grid_direction.size;
#endif

      //#pragma omp parallel default(none) firstprivate(count_directions) shared(sorted_graph, sorted_id_bound_face, inner_bound_code, vec_x0, grid, norm)
      {
        const IdType count_cells = grid.size;

        Type loc_norm = -1;
        std::vector<Type> *inter_coef = &grid.inter_coef_all[omp_get_thread_num()]; ///< указатель на коэффициенты интерполяции по локальному для потока направлению

        //#pragma omp for
        for (IdType num_direction = 0; num_direction < count_directions; ++num_direction) {

          /*---------------------------------- FOR по граничным граням----------------------------------*/

          for (auto num_face : sorted_id_bound_face[num_direction]) {
            const IdType neigh_id = grid.faces[num_face].geo.id_r;                    //т.к. мы проходим границу здесь будет граничный признак
            (*inter_coef)[num_face] = spec::BoundaryConditions(neigh_id, frq0, frq1); //значение на грани ( или коэффициенты интерполяции)
          }

          // индексация по массиву определяющих гранях (конвеерная т.к. заранее не известны позиции точек)
          const Type *s = vec_x0[num_direction].s.data();
          const face_loc_id_t *in_face = vec_x0[num_direction].in_face_id.data();
          /*---------------------------------- далее FOR по неосвященным граням----------------------------------*/
          for (auto fc_pair : sorted_graph[num_direction]) {

            const IntId num_cell = fc_pair.cell;
            const uint32_t num_loc_face = fc_pair.loc_face;
            elem_t *cell = &grid.cells[num_cell];
            const Vector3 &x = cell->geo.center; // vec_x[num_cell].x[num_loc_face][num_node];
#if 0
          const Type S = grid.scattering[num_direction * count_cells + num_cell];
          Type k;
          const Type rhs = GetRhs(x, S, *cell, k);
          const face_loc_id_t id_in_faces = *(in_face);
          ++in_face;
          alignas(32) Type I0[4] =
              {(*inter_coef)[cell->geo.id_faces[id_in_faces.a]],
               (*inter_coef)[cell->geo.id_faces[id_in_faces.b]],
               (*inter_coef)[cell->geo.id_faces[id_in_faces.c]],
               0};
          (*inter_coef)[cell->geo.id_faces[num_loc_face]] = GetIllum(I0, s, k, rhs);
          s += NODE_SIZE;
#else
            const Type S = grid.scattering[num_direction * count_cells + num_cell];
            const face_loc_id_t id_in_faces = *(in_face);
            ++in_face;
            alignas(32) Type I0[4] =
                {(*inter_coef)[cell->geo.id_faces[id_in_faces.a]],
                 (*inter_coef)[cell->geo.id_faces[id_in_faces.b]],
                 (*inter_coef)[cell->geo.id_faces[id_in_faces.c]],
                 0};
            Type I = 0;
            // структура аналогичная  ::trace::GetLocNodes(...)
            for (ShortId num_node = 0; num_node < NODE_SIZE; ++num_node) {
              const Vector3 &dir = grid_direction.directions[num_direction].dir;
              I += spec::GetIllum(dir, x, s[num_node], I0[num_node], S, frq0, frq1, *cell);
            }
            // WRITE_LOG("fI=%lf, \n", I);
            s += NODE_SIZE;
            //записываем значение коэффициентов на грани
            (*inter_coef)[cell->geo.id_faces[num_loc_face]] = I / NODE_SIZE;
#endif
          }

          /*---------------------------------- конец FOR по ячейкам----------------------------------*/
          loc_norm = spec::ReCalcIllum(num_direction, num_frq, *inter_coef, grid);
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

      // WRITE_LOG("End iter number#: %d, norm=%.16lf, time= %lf\n", iter, norm,
      //           (double)tick::duration_cast<tick::milliseconds>(tick::steady_clock::now() - start_clock).count() / 1000.);
      iter++;
    } while (norm > _solve_mode.accuracy && iter < _solve_mode.max_number_of_iter);

    res_spectrum[num_frq] = get_full_illum(grid_direction.size, num_frq, grid);

    // res_spectrum[num_frq] = 0;
    // for (size_t j = 0; j < CELL_SIZE; j++) {
    //   constexpr int cellid = 4719;
    //   res_spectrum[num_frq] += grid.Illum[CELL_SIZE * (0 + cellid) + j]; // get_full_illum(grid_direction.size, num_frq, grid);
    // }
    // res_spectrum[num_frq] /= CELL_SIZE;

  } // num_freq

  std::vector<Type> frqs(spectrum.size() - 1);
  for (size_t i = 0; i < frqs.size(); i++) {
    frqs[i] = (spectrum[i + 1] + spectrum[i]) / 2.0;
  }

  WriteRes(glb_files.solve_address + "SumIllum.txt", frqs, res_spectrum);

  std::vector<Type> source(spectrum.size() - 1);
  for (size_t i = 0; i < source.size(); i++) {
    source[i] = exp(B_Plank_log(5000.0, spectrum[i + 1], spectrum[i]));
  }
  WriteRes(glb_files.solve_address + "source.txt", source);

  return e_completion_success;
}

#endif