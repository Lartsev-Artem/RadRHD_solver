#if defined SOLVERS && defined ILLUM
#include "illum_params.h"

#include "illum_integrate.h"
#include "illum_utils.h"

void illum::GetEnergy(const grid_directions_t &grid_direction, grid_t &grid) {

#pragma omp parallel default(none) shared(grid_direction, grid)
  {
#pragma omp for
    for (int i = 0; i < grid.size; i++) {
      elem_t &el = grid.cells[i];
      el.illum_val.energy = direction_integrator::IntegrateByCell(el.illum_val.illum, grid_direction);
    }
  }
}

void illum::GetStream(const grid_directions_t &grid_direction, grid_t &grid) {

#pragma omp parallel default(none) shared(grid_direction, grid)
  {
#pragma omp for
    for (int i = 0; i < grid.cells.size(); i++) {
      elem_t &el = grid.cells[i];

      Vector3 Stream[CELL_SIZE];
      direction_integrator::IntegrateByFace3(el.illum_val.illum, grid_direction, Stream);

      el.illum_val.stream = GetAverageByCell(Stream); //поток

      el.illum_val.div_stream = 0;
      for (int j = 0; j < CELL_SIZE; j++) {
        geo_face_t *geo_f = &grid.faces[el.geo.id_faces[j]].geo;
        if (el.geo.sign_n[j]) {
          el.illum_val.div_stream += Stream[j].dot(geo_f->n) * geo_f->S;
        } else {
          el.illum_val.div_stream -= Stream[j].dot(geo_f->n) * geo_f->S;
        }
      }
      el.illum_val.div_stream /= el.geo.V;
    }
  }
}

void illum::GetImpuls(const grid_directions_t &grid_direction, grid_t &grid) {

#pragma omp parallel default(none) shared(grid_direction, grid)
  {
#pragma omp for
    for (int i = 0; i < grid.cells.size(); i++) {
      elem_t &el = grid.cells[i];

      Matrix3 Impuls[CELL_SIZE];
      direction_integrator::IntegrateByFace9(el.illum_val.illum, grid_direction, Impuls);

      el.illum_val.impuls = GetAverageByCell(Impuls); //импульс

      el.illum_val.div_impuls = Vector3::Zero();
      for (int j = 0; j < CELL_SIZE; j++) {
        geo_face_t *geo_f = &grid.faces[el.geo.id_faces[j]].geo;
        if (el.geo.sign_n[j]) {
          el.illum_val.div_impuls += Impuls[j] * (geo_f->n) * geo_f->S;
        } else {
          el.illum_val.div_impuls -= Impuls[j] * (geo_f->n) * geo_f->S;
        }
      }
      el.illum_val.div_impuls /= el.geo.V;
    }
  }
}

#endif //! defined SOLVERS && defined ILLUM