#if 0
int TestDivStream(const std::vector<Vector3> &centers_face, grid_t &grid) {
#pragma omp parallel default(none) shared(centers_face, grid)
  {
#pragma omp for
    for (int i = 0; i < grid.cells.size(); i++) {
      elem_t &el = grid.cells[i];

      Vector3 Stream[CELL_SIZE];
      for (int j = 0; j < CELL_SIZE; j++) {
        Vector3 x = centers_face[i * CELL_SIZE + j];
        Stream[j] = Vector3(3 * x[0], 0, 0);
      }

      GET_FACE_TO_CELL(el.illum_val.stream, Stream, Vector3::Zero());

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
  return 0;
}
#endif