
#include "trace_to_face.h"

void trace::GetBoundGraphFaces(const grid_t &grid, const std::vector<IntId> &neighbours,
                               const std::vector<bits_flag_t> &face_states, const std::vector<IntId> &sorted_id_cell,
                               std::vector<IntId> &graph_bound_faces) {
  graph_bound_faces.clear();
  graph_bound_faces.reserve(grid.size / 10);

  for (auto num_cell : sorted_id_cell) {
    for (ShortId num_out_face = 0; num_out_face < CELL_SIZE; ++num_out_face) {

      // если эта грань входящая и граничная
      if (CHECK_BIT(face_states[num_cell], num_out_face) == e_face_type_in) {
        if (neighbours[num_cell * CELL_SIZE + num_out_face] < 0) ///< сосед к текущей грани
        {
          graph_bound_faces.push_back(grid.cells[num_cell].geo.id_faces[num_out_face]); //номера граней на сетке
        }
      }
    }
  }
  // files_sys::bin::WriteSimple(glb_files.illum_geo_address + "faces_bound" + std::to_string(dir) + ".bin", graph_bound_faces);
  return;
}

void trace::GetGraphFaces(const grid_t &grid,
                          const std::vector<bits_flag_t> &face_states, const std::vector<IntId> &sorted_id_cell,
                          std::vector<IntId> &graph_cell_faces) {

  graph_cell_faces.clear();
  graph_cell_faces.reserve(2 * grid.faces.size());

  for (auto num_cell : sorted_id_cell) {
    // расчитываем излучения на выходящих гранях
    for (ShortId num_out_face = 0; num_out_face < CELL_SIZE; ++num_out_face) {
      // если эта грань входящая и граничная, то пропускаем её
      if (CHECK_BIT(face_states[num_cell], num_out_face) == e_face_type_in) {
        continue;
      }
      ///\note: в теории можно хранить ячейку и номер в локальной нумерации, а ссылку получать в расчёте
      graph_cell_faces.push_back(num_cell);
      graph_cell_faces.push_back(grid.cells[num_cell].geo.id_faces[num_out_face]); //это вместо graph
    }                                                                              // num_out_face
  }
  /*---------------------------------- конец FOR по ячейкам----------------------------------*/

  // files_sys::bin::WriteSimple(glb_files.illum_geo_address + "faces_in_body" + std::to_string(num_direction) + ".bin", faces_in_body);
  // files_sys::bin::WriteSimple(glb_files.illum_geo_address + "faces_body" + std::to_string(num_direction) + ".bin", faces_body);

  return;
}