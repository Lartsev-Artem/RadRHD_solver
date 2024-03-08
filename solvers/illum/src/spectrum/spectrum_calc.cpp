#if defined SOLVERS && defined ILLUM && SPECTRUM
#include "spectrum_calc.h"
#if defined TRANSFER_CELL_TO_FACE && defined SEPARATE_GPU

#include "illum_calc_gpu_async.h"
#include "spectrum_utils.h"

#include "cuda_interface.h"

#include "global_consts.h"
#include "illum_utils.h"
#include "plunk.h"

static constexpr int projection_direction = 66;
bool log_enable = 0;

int illum::spectrum::CalculateSpectrum(const grid_directions_t &grid_direction,
                                       const std::vector<std::vector<IntId>> &inner_bound_code,
                                       const std::vector<align_cell_local> &vec_x0,
                                       const std::vector<std::vector<graph_pair_t>> &sorted_graph,
                                       const std::vector<std::vector<IntId>> &sorted_id_bound_face,
                                       grid_t &grid) {

  const int count_frequencies = grid.frq_grid.size() - 1; ///< число разбиений по частоте
  grid.spectrum.assign(grid.spectrum.size(), 0);

  for (int num_frq = 0; num_frq < count_frequencies; num_frq++) {
    // for (int num_frq = 0; num_frq < 1; num_frq++) {

    grid.cur_frq_id = num_frq;
    const Type frq0 = grid.frq_grid[grid.cur_frq_id];
    const Type frq1 = grid.frq_grid[grid.cur_frq_id + 1];
    Type val = exp(B_Plank_log(40000000, frq1, frq0) - LOG(kRadiation));
    log_spectrum("frq[%d]= %lf %lf, Bnd=%lf\n", grid.cur_frq_id, frq0, frq1, val);
    if (val < 1e-200) {
      log_spectrum("Skip frq %d\n", num_frq);
      continue;
    }
    set_boundary_value(val);

    illum::separate_gpu::CalculateIllum(grid_direction, inner_bound_code,
                                        vec_x0, sorted_graph, sorted_id_bound_face, grid);

    cuda::interface::CudaWait();
    cuda::interface::CudaSyncStream(cuda::e_cuda_scattering_1);
    cuda::interface::CudaSyncStream(cuda::e_cuda_params);
    grid.spectrum[num_frq] = get_full_illum(projection_direction, grid);
    WRITE_LOG("end frq %d: %.16lf\n", num_frq, grid.spectrum[num_frq]);
  }

  log_spectrum("dir: %lf %lf %lf\n", grid_direction.directions[projection_direction].dir[0],
               grid_direction.directions[projection_direction].dir[1],
               grid_direction.directions[projection_direction].dir[2]);

  return 0;
}
#endif
#endif