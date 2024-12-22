#include "convert_bin_to_geo.h"

void ConvertBinToGeo(std::vector<IntId> &neighbours_id_faces,
                     std::vector<Normals> &normals, std::vector<Type> &areas_faces,
                     std::vector<Type> &volume, std::vector<Vector3> &centers,
                     std::vector<face_t> &faces, std::vector<elem_t> &cells) {

  const int N = centers.size();
  cells.resize(N);

  int cc = 0;
  for (int i = 0; i < N * CELL_SIZE; i++) {
    int idx = neighbours_id_faces[i];

    if (idx != -10) //эта грань еще не обработана
    {
      face_t f;
      f.geo.id_l = i / CELL_SIZE; //ячейка

      f.geo.n = normals[i / CELL_SIZE].n[i % CELL_SIZE];
      f.geo.S = areas_faces[i];

      cells[i / CELL_SIZE].geo.sign_n.set_sign(i % CELL_SIZE, true);
      cells[i / CELL_SIZE].geo.id_faces[i % CELL_SIZE] = cc;

      neighbours_id_faces[i] = -10;
      if (idx >= 0) {
        f.geo.id_r = idx / CELL_SIZE; //сосед

        cells[idx / CELL_SIZE].geo.sign_n.set_sign(idx % CELL_SIZE, false);
        cells[idx / CELL_SIZE].geo.id_faces[idx % CELL_SIZE] = cc;

        neighbours_id_faces[idx] = -10;
      } else {
        f.geo.id_r = idx;                    // код границы
        cells[i / CELL_SIZE].geo.is_bound++; //сразу определяет число граничных граней (для углов)
      }

      faces.push_back(f); // как потом искать с ячейками?
      cc++;
    }
  }

  for (int i = 0; i < N; i++) {
    cells[i].geo.center = centers[i];
    cells[i].geo.V = volume[i];
  }
}