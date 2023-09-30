#if defined SOLVERS && defined ILLUM
#include "scattering.h"

#include "illum_utils.h"

Type illum::scattering::GetIntScatteringByDirection(const int num_cell, const Vector3 &direction, const int size_one_direction, const Type *illum_old,
                                                    const grid_directions_t &grid_direction) {

  // const int size_one_direction = illum_old.size() / N_dir;
  Type S = 0;
  for (int num_direction = 0; num_direction < grid_direction.size; num_direction++) {

    Type I = GetAverageByCell(illum_old + (num_direction * size_one_direction + num_cell));
    S += Gamma(grid_direction.directions[num_direction].dir, direction) * I * grid_direction.directions[num_direction].area;
  }

  return S / grid_direction.full_area; // было *4PI, но из-за нормировки Gamma разделили на 4PI
}

void illum::scattering::CalculateIntCPU(const grid_directions_t &grid_direction, grid_t grid) {

#pragma omp parallel default(none) shared(grid, grid_direction)
  {
    Vector3 direction;
    const int N = grid.size;
#pragma omp for
    for (int num_direction = 0; num_direction < grid_direction.size; ++num_direction) {

      direction = grid_direction.directions[num_direction].dir;
      for (int cell = 0; cell < N; cell++) {
        grid.scattering[num_direction * N + cell] = GetIntScatteringByDirection(CELL_SIZE * cell, direction, N * CELL_SIZE, grid.Illum, grid_direction);
      }
    }
  }
}

#endif //! defined SOLVERS && defined ILLUM