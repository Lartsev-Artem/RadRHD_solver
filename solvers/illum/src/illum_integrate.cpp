#if defined SOLVERS && defined ILLUM
#include "illum_integrate.h"
#include "illum_utils.h"

namespace integ = illum::direction_integrator;

Type integ::IntegrateByCell(const std::vector<Type> &Illum, const grid_directions_t &grid_direction) {

  Type res = 0;
  int i = 0;
  for (auto &dir : grid_direction.directions) {
    res += GetAverageByCell(&Illum[i * CELL_SIZE]) * dir.area;
    i++;
  }

  return res / grid_direction.full_area;
}

void integ::IntegrateByFace3(const std::vector<Type> &Illum, const grid_directions_t &grid_direction, Vector3 *stream_face) {
  int i = 0;
  for (int f = 0; f < CELL_SIZE; f++) {
    stream_face[f] = Vector3::Zero();
  }

  for (auto &dir : grid_direction.directions) {
    for (int f = 0; f < CELL_SIZE; f++) {
      stream_face[f] += Illum[i++] * dir.area * dir.dir;
    }
  }

  for (int f = 0; f < CELL_SIZE; f++) {
    stream_face[f] /= grid_direction.full_area;
  }
}

void integ::IntegrateByFace9(const std::vector<Type> &Illum, const grid_directions_t &grid_direction, Matrix3 *impuls_face) {

  int i = 0;
  for (int f = 0; f < CELL_SIZE; f++) {
    impuls_face[f] = Matrix3::Zero();
  }

  for (auto &dir : grid_direction.directions) {
    for (int f = 0; f < CELL_SIZE; f++) {

      for (int h = 0; h < 3; h++)
        for (int k = 0; k < 3; k++) {
          impuls_face[f](h, k) += dir.dir[h] * dir.dir[k] * (Illum[i] * dir.area);
        }
      i++;
    }
  }

  for (int f = 0; f < CELL_SIZE; f++) {
    impuls_face[f] /= grid_direction.full_area;
  }
}

#endif //! defined SOLVERS && defined ILLUM