#if 0 //здесь есть схема расчета полного излучения. С хранением всего спектра
// => верного расчёта интеграла рассеяния и численных интегралов по частоте
// однако реализация требует N*M*F памяти на узел, что допускает только 
// грубые сетки + нет окончательного модуля к видеокарте. 
// Схема запускалась только без рассеяния и требует отладки

/// \todo: нужна проверка отправлений. Тестовый запуск без видеокарты. Все таки подумать над Illum_loc
#include "solvers_config.h"
#if defined SPECTRUM && defined TRANSFER_CELL_TO_FACE
#include "spec_all.h"

#ifdef USE_CUDA
#include "cuda_interface.h"
#endif
#include <omp.h>

#define USE_ONE_DIRECTION

#include "global_value.h"
#include "plunk.h"

#ifdef USE_TRACE_THROUGH_INNER_BOUNDARY
const IntId *code_bound = inner_bound_code[num_direction].data();
#pragma error "need recalc code bound from cell to face"
#endif

Type illum::spec::ReCalcIllum(const IdType num_dir, const std::vector<std::vector<Type>> &inter_coef, grid_t &grid, const IdType dir_disp) {
  Type norm = -1;
  const IdType shift_dir = (grid.size_frq * (num_dir + dir_disp));
  const IdType shift_cell = grid.size_frq * grid.size_dir;

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

  return norm;
}

#include "writer_txt.h"
static void WriteRes(const std::string &name_file, std::vector<Type> &y) {
  files_sys::txt::WriteSimple(name_file, y);
}

static void WriteRes(const std::string &name_file, std::vector<Type> &x, std::vector<Type> &y) {
  files_sys::txt::WriteSimple(name_file, x, y);
}

#ifndef SEPARATE_GPU
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

#else // SEPARATE_GPU
namespace tick = std::chrono;
#include "cuda_multi_interface.h"
#include "illum_mpi_sender.h"
namespace cuda_sep = cuda::interface::separate_device;

#include "plunk.h"
static Type GetRhs(const Type S, Type frq1, Type frq0, elem_t &cell, Type &k) {
  full_phys_data_t *phys = &cell.cell_data;

  Type betta;
  if (LIKELY(phys->vel > kC_LightInv)) {
    betta = (get_scat_coef(0.5 * (frq1 + frq0), phys->vel, phys->cosf, phys->lorenz) / (kM_hydrogen * kDist)) * phys->val->d;
  } else {
    betta = (get_scat_coef(0.5 * (frq1 + frq0)) / (kM_hydrogen * kDist)) * phys->val->d;
  }

  Type alpha = phys->alpha;
  Type Q = B_Plank(phys->T, phys->logT, frq1, frq0); // Q=alpha*Ie

  k = alpha + betta;
  return (alpha * Q + betta * S) / k;
}
#pragma GCC target("avx2")
#pragma GCC optimize("O3")

#include <bits/stdc++.h>
#include <x86intrin.h>
static inline __m256d _mm256_exp_pd(__m256d x) {
  alignas(32) Type X[4];
  _mm256_store_pd(X, x);
  return _mm256_setr_pd(
      exp(X[0]),
      exp(X[1]),
      exp(X[2]),
      0);
}
static Type GetIllum2(const Type *I0, const Type *s, const Type k, const Type rhs) { //это ок

  alignas(32) Type Icur[4];
  __m256d RHS = _mm256_set1_pd(rhs);
  __m256d S = _mm256_load_pd(s);
  __m256d I_0 = _mm256_load_pd(I0);

  __m256d EXP = _mm256_exp_pd(_mm256_mul_pd(_mm256_set1_pd(-k), S));
  __m256d I = _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(I_0, RHS), EXP), RHS);
  _mm256_store_pd(Icur, _mm256_max_pd(_mm256_set1_pd(0), I));
  return (Icur[0] + Icur[1] + Icur[2]) / 3.;
}
#pragma GCC pop("O3");

int illum::spec::CalculateIllumFaceMpi(const grid_directions_t &grid_direction,
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

  for (auto &el : grid.cells) {
    el.cell_data.Init(&el.phys_val);
  }

  do {
    auto start_clock = tick::steady_clock::now();
    norm = -1;

    ///\note здесь ждем конкретную частоту
    // if (section_1.requests_send[0] != MPI_REQUEST_NULL) {
    //   MPI_Waitall(section_1.requests_send.size(), section_1.requests_send.data(), MPI_STATUSES_IGNORE);
    // }
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

#pragma omp parallel default(none) firstprivate(count_directions, myid, np, local_disp, local_size, count_frequencies) \
    shared(grid_direction, sorted_graph, sorted_id_bound_face, inner_bound_code, vec_x0, grid, norm, section_1)
    {

      const int count_th = omp_get_num_threads();
      const int num_th = omp_get_thread_num();

      const IdType count_cells = grid.size;

      Type loc_norm = -1;
      std::vector<std::vector<Type>> *inter_coef = &grid.inter_coef_all[omp_get_thread_num()]; ///< указатель на коэффициенты интерполяции по локальному для потока направлению

#pragma omp single // nowait
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

          const Type *frq = grid.spectrum.data() - 1; //-1 т.к. на первой итерации сдвиг
          std::vector<Type> *init_frq_coef = &(*inter_coef)[num_face];
          for (IdType num_frq = 0; num_frq < count_frequencies - 1; num_frq++) {
            Type frq0 = *(++frq);
            /// \todo можно прокинуть общий neigh_id и вызывать частотный инициализатор без доп ветвления
            (*init_frq_coef)[num_frq] = illum::spec::BoundaryConditions(neigh_id, frq0, *frq); //значение на грани ( или коэффициенты интерполяции)
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
          cell->cell_data.InitDirection(grid_direction.directions[num_direction].dir);
          std::vector<Type> *frq_coef_A = &(*inter_coef)[cell->geo.id_faces[id_in_faces.a]];
          std::vector<Type> *frq_coef_B = &(*inter_coef)[cell->geo.id_faces[id_in_faces.b]];
          std::vector<Type> *frq_coef_C = &(*inter_coef)[cell->geo.id_faces[id_in_faces.c]];
          std::vector<Type> *frq_coef_Res = &(*inter_coef)[cell->geo.id_faces[num_loc_face]];
          const Type *frq = grid.spectrum.data() - 1; //-1 т.к. на первой итерации сдвиг

          for (IdType num_frq = 0; num_frq < count_frequencies - 1; num_frq++) {
            Type frq0 = *(++frq);
            Type k;
            const Type rhs = GetRhs(S[num_frq], frq0, *frq, *cell, k);

            alignas(32) Type I0[4] = {(*frq_coef_A)[num_frq], (*frq_coef_B)[num_frq], (*frq_coef_C)[num_frq], 0};
            (*frq_coef_Res)[num_frq] = GetIllum2(I0, s, k, rhs);
          }
          /*---------------------------------- конец FOR по частотам----------------------------------*/
          s += NODE_SIZE;
        }
        /*---------------------------------- конец FOR по ячейкам----------------------------------*/
        loc_norm = spec::ReCalcIllum(num_direction, *inter_coef, grid, local_disp);

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
      cuda_sep::CalculateIntScatteringAsync(grid_direction, grid, 0, local_size, cuda::e_cuda_scattering_1);
    }

    MPI_Wait(&rq_norm, MPI_STATUS_IGNORE);
    norm = fabs(glob_norm);

    WRITE_LOG("End iter number#: %d, norm=%.16lf, time= %lf\n", iter, norm,
              (double)tick::duration_cast<tick::milliseconds>(tick::steady_clock::now() - start_clock).count() / 1000.);
    iter++;
  } while (norm > _solve_mode.accuracy && iter < _solve_mode.max_number_of_iter);

  return e_completion_success;
}
#endif
#endif

#endif