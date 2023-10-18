/**
 * @file json_struct.cpp
 * @brief Пользовательские методы для json структур
 */

#include "json_struct.h"

#include "global_value.h"

global_files_t glb_files; ///< глобальная структура с путями к данным

global_files_t::global_files_t(const std::string &BASE_ADDRESS, const std::string &address_illum_geo_file) {

  base_address = BASE_ADDRESS;
  illum_geo_address = address_illum_geo_file;
  graph_address = F_GRAPH;
  solve_address = F_SOLVE;
  hllc_init_value = F_INIT_HLLC;
  Build();
}

void global_files_t::Build() {
  // name_file_hllc_set = base_address + F_HLLC_SET;

  //---------Файлы расчётных данных(отдельные файлы по направлениям)---------//
  name_file_state_face = illum_geo_address + F_STATE_FACE;
  name_file_x0_loc = illum_geo_address + F_X0_LOC;
  name_file_x = illum_geo_address + F_X;

  //--------Файлы трассировки сквозь внутреннюю границу------//
  name_file_dist_try = base_address + F_DIST_TRY;
  name_file_id_try = base_address + F_ID_TRY;
  name_file_res = graph_address + F_RES;
  name_file_neigh = base_address + F_NEIGHBOR;

#if 0
		//--------Файлы сдвигов по направлениям в расчётных файлах------------------//
		const std::string name_file_shift_out = glb_files.base_address + "ShiftOut";
		const std::string name_file_shift_res = glb_files.base_address + "ShiftRes";
		const std::string name_file_shift_x0 = glb_files.base_address + "ShiftX0";
		const std::string name_file_shift_try = glb_files.base_address + "ShiftTry";
#endif

  //-------------------Файлы геометрии-------------------//
  name_file_geometry_faces = base_address + F_GEO_FACES;
  name_file_geometry_cells = base_address + F_GEO_CELLS;
}

#ifdef DEBUG
void global_files_t::print() {
  std::string *str = (std::string *)this;
  while (str < (std::string *)this + sizeof(global_files_t) / sizeof(std::string)) {
    std::cout << *str++ << '\n';
  }
}
#endif
