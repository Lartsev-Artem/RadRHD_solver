#include "reader_bin.h"

#include <iostream>

#include "global_def.h"
#include "global_value.h"

#include "mpi_shifts.h"

int files_sys::bin::ReadNormals(const std::string &name_file_normals, std::vector<Normals> &normals) {
  FILE *f;
  OPEN_FILE(f, name_file_normals.c_str(), "rb");

  int n;
  if (fread(&n, sizeof(int), 1, f) != 1) {
    return e_completion_fail;
  }
  normals.resize(n);

  Normals norm;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < CELL_SIZE; j++)
      if (fread(&norm.n[j], sizeof(Vector3), 1, f) != 1) {
        return e_completion_fail;
      }

    normals[i] = norm;
  }
  fclose(f);

  return e_completion_success;
}

int files_sys::bin::ReadData(const size_t class_file_vtk, const std::string &main_dir,
                             std::vector<Type> &density, std::vector<Type> &absorp_coef,
                             std::vector<Type> &rad_en_loose_rate,
                             std::vector<Vector3> &velocity, std::vector<Type> &pressure,
                             const bool is_print /*=false*/) {
  switch (class_file_vtk) {
  case e_grid_cfg_radiation:
    ReadSimple(main_dir + F_DENSITY, density);
    ReadSimple(main_dir + F_ABSORPCOEF, absorp_coef);
    ReadSimple(main_dir + F_RADLOOSERATE, rad_en_loose_rate);
    break;

  case e_grid_cfg_full_init:
    ReadSimple(main_dir + F_DENSITY, density);
    ReadSimple(main_dir + F_ABSORPCOEF, absorp_coef);
    ReadSimple(main_dir + F_RADLOOSERATE, rad_en_loose_rate);

    ReadSimple(main_dir + F_VELOCITY, velocity);
    ReadSimple(main_dir + F_PRESSURE, pressure);
    break;

  default:
    printf("Grid without data\n");
    return e_completion_success;
  }

  if (is_print) {
    std::cout << "density_Size: " << density.size() << '\n';
    std::cout << "absorp_coef_Size: " << absorp_coef.size() << '\n';
    std::cout << "Q_Size: " << rad_en_loose_rate.size() << '\n';
  }

  const size_t a = density.size();
  const size_t b = rad_en_loose_rate.size();
  const size_t c = absorp_coef.size();
  if (a == b && a == c && b == c) {
    return e_completion_success;
  }

  RETURN_ERR("Error size data\n");
}

#if 0 // def ILLUM
#define READ_FILE(name_file, data, value)       \
  {                                             \
    FILE *f;                                    \
    OPEN_FILE(f, name_file., "rb");             \
    int n;                                      \
    fread(&n, sizeof(int), 1, f);               \
    data.resize(n);                             \
    for (auto &el : data) {                     \
      fread(&el.value, sizeof(el.value), 1, f); \
    }                                           \
    fclose(f);                                  \
  }

//! Не использовать такую структура чтения!. resize cells не нужен!!!. При
//! излучении данные распределены по массивам?
int files_sys::bin::ReadData(const solve_mode_t &mode, const std::string &main_dir,
                             grid_t &grid) {
  switch (mode.class_vtk) {
  case 1: // test grid
    READ_FILE((main_dir + "alpha.bin").c_str(), grid.cells, phys_val.d);
    READ_FILE((main_dir + "alpha.bin").c_str(), grid.cells,
              illum_val.absorp_coef);
    READ_FILE((main_dir + "Q.bin").c_str(), grid.cells,
              illum_val.rad_en_loose_rate);
    break;

  case 2: // main grid to illum
    READ_FILE((main_dir + "density.bin").c_str(), grid.cells, phys_val.d);
    READ_FILE((main_dir + "AbsorpCoef.bin").c_str(), grid.cells,
              illum_val.absorp_coef);
    READ_FILE((main_dir + "radEnLooseRate.bin").c_str(), grid.cells,
              illum_val.rad_en_loose_rate);
    break;

  case 3: // full_test_grid
    READ_FILE((main_dir + "density.bin").c_str(), grid.cells, phys_val.d);
    READ_FILE((main_dir + "alpha.bin").c_str(), grid.cells,
              illum_val.absorp_coef);
    READ_FILE((main_dir + "Q.bin").c_str(), grid.cells,
              illum_val.rad_en_loose_rate);
    READ_FILE((main_dir + "velocity.bin").c_str(), grid.cells, phys_val.v);
    READ_FILE((main_dir + "pressure.bin").c_str(), grid.cells, phys_val.p);
    break;
  default:
    printf("Grid without data\n");
    return 0;
  }

  return 0;
}
#endif
#ifdef ILLUM
int files_sys::bin::ReadRadiationTrace(const int count_dir, const global_files_t &gbl_files,
                                       std::vector<BasePointTetra> &vec_x,
                                       std::vector<std::vector<State>> &face_states,
                                       std::vector<std::vector<cell_local>> &vec_x0,
                                       std::vector<std::vector<IntId>> &sorted_id_cell,
                                       std::vector<std::vector<IntId>> &inner_bound_code) {
  if (ReadSimple(gbl_files.name_file_x, vec_x))
    return e_completion_fail;

  std::vector<IdType> disp;
  std::vector<IdType> send;

  int np = get_mpi_np();
  int myid = get_mpi_id();

  GetDisp(np, count_dir, disp);
  GetSend(np, count_dir, send);

  vec_x0.resize(send[myid]);

  sorted_id_cell.resize(send[myid]);
  face_states.resize(send[myid]);
  inner_bound_code.resize(send[myid]);

  for (int i = 0; i < send[myid]; i++) {

    if (ReadSimple(gbl_files.name_file_x0_loc + std::to_string(disp[myid] + i) + ".bin", vec_x0[i]))
      return e_completion_fail;

    if (ReadSimple(gbl_files.graph_address + F_GRAPH + std::to_string(disp[myid] + i) + ".bin", sorted_id_cell[i]))
      return e_completion_fail;

    if (ReadSimple(gbl_files.name_file_state_face + std::to_string(disp[myid] + i) + ".bin", face_states[i]))
      return e_completion_fail;
  }

#ifdef USE_TRACE_THROUGH_INNER_BOUNDARY
  for (int i = 0; i < send[myid]; i++) {
    if (ReadSimple(gbl_files.name_file_res + std::to_string(disp[myid] + i) + ".bin", inner_bound_code[i])) {
      WRITE_LOG("WARNING!!! inner bound[%d] didn't read\n", i);
      inner_bound_code.clear();
      break;
    }
  }
#endif

  return e_completion_success;
}

int files_sys::bin::ReadRadiationFaceTrace(const int count_dir, const global_files_t &gbl_files,
                                           std::vector<BasePointTetra> &vec_x,
                                           std::vector<std::vector<cell_local>> &vec_x0,
                                           std::vector<std::vector<graph_pair_t>> &sorted_graph,
                                           std::vector<std::vector<IntId>> &sorted_id_bound_face,
                                           std::vector<std::vector<IntId>> &inner_bound_code) {
  if (ReadSimple(gbl_files.name_file_x, vec_x))
    return e_completion_fail;

  std::vector<IdType> disp;
  std::vector<IdType> send;

  int np = get_mpi_np();
  int myid = get_mpi_id();

  GetDisp(np, count_dir, disp);
  GetSend(np, count_dir, send);

  vec_x0.resize(send[myid]);

  sorted_id_bound_face.resize(send[myid]);
  sorted_graph.resize(send[myid]);
  inner_bound_code.resize(send[myid]);

  for (int i = 0; i < send[myid]; i++) {

    if (ReadSimple(gbl_files.name_file_x0_loc + std::to_string(disp[myid] + i) + ".bin", vec_x0[i]))
      return e_completion_fail;

    if (ReadSimple(gbl_files.graph_address + F_GRAPH_BOUND_FACE + std::to_string(disp[myid] + i) + ".bin", sorted_id_bound_face[i]))
      return e_completion_fail;

    if (ReadSimple(gbl_files.graph_address + F_GRAPH_BODY_FACE + std::to_string(disp[myid] + i) + ".bin", sorted_graph[i]))
      return e_completion_fail;
  }

#ifdef USE_TRACE_THROUGH_INNER_BOUNDARY
  for (int i = 0; i < send[myid]; i++) {
    if (ReadSimple(gbl_files.name_file_res + std::to_string(disp[myid] + i) + ".bin", inner_bound_code[i])) {
      WRITE_LOG("WARNING!!! inner bound[%d] didn't read\n", i);
      inner_bound_code.clear();
      break;
    }
  }
#endif

  return e_completion_success;
}
#endif //! ILLUM
