#if defined SOLVERS && defined ILLUM && SPECTRUM
#include "spectrum_calc.h"
#if defined TRANSFER_CELL_TO_FACE && defined SEPARATE_GPU

#include "illum_calc_gpu_async.h"
#include "spectrum_utils.h"

#include "cuda_interface.h"

int illum::spectrum::CalculateSpectrum(const grid_directions_t &grid_direction,
                                       const std::vector<std::vector<IntId>> &inner_bound_code,
                                       const std::vector<align_cell_local> &vec_x0,
                                       const std::vector<std::vector<graph_pair_t>> &sorted_graph,
                                       const std::vector<std::vector<IntId>> &sorted_id_bound_face,
                                       grid_t &grid) {

  const int count_frequencies = grid.frq_grid.size() - 1; ///< число разбиений по частоте

  for (int num_frq = 0; num_frq < count_frequencies; num_frq++) {

    grid.cur_frq_id = num_frq;

    illum::separate_gpu::CalculateIllum(grid_direction, inner_bound_code,
                                        vec_x0, sorted_graph, sorted_id_bound_face, grid);

    cuda::interface::CudaWait();
    cuda::interface::CudaSyncStream(cuda::e_cuda_scattering_1);
    cuda::interface::CudaSyncStream(cuda::e_cuda_params);
    grid.spectrum[num_frq] = get_full_illum(0, grid);
  }

  return 0;
}
#endif
#endif