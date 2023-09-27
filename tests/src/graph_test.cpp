/**
 * @file graph_test.cpp
 * @brief Тест для модуля построителя графов
 * @version 0.1
 * @date 2023-09-26
 *
 * \details тест из файла vtk и сферы направлений строит граф обхода и создает новую сетку
 * со скалярными полями где значение --- позиция ячейки в графе
 *
 * \note требуется поддержка VTK
 *
 */

#include "build_internal_format.h"
#include "graph_main.h"

#include "global_value.h"
#include "reader_bin.h"
#include "writer_bin.h"

#include "file_module.h"

#ifdef USE_VTK
int main(int argc, char **argv) {

  std::string prj_dir = "/home/artem/projects/solver/";

  glb_files.base_address = prj_dir + " tests/build/";
  glb_files.graph_address = prj_dir + "tests/build/graph/";
  glb_files.name_file_sphere_direction = prj_dir + "tests/data/surface_26_dir.txt";
  glb_files.name_file_vtk = prj_dir + "tests/data/sphere.vtk";
  glb_files.Build();
  std::string file_vtk_out = glb_files.base_address + "graph_test.vtk";

  fs::create_directory(glb_files.base_address);
  fs::create_directory(glb_files.graph_address);
  fs::copy_file(glb_files.name_file_vtk, file_vtk_out, fs::copy_options::overwrite_existing);

  if (BuildDataFromVTK(glb_files) != e_completion_success) {
    return e_completion_fail;
  }

  if (graph::RunGraphModule() != e_completion_success) {
    return e_completion_fail;
  }

  std::vector<int> graph;
  std::vector<int> graph_order;

  for (const auto &entry : fs::directory_iterator(glb_files.graph_address)) //по всем файлам в директории
  {
    if (entry.is_regular_file()) {
      files_sys::bin::ReadSimple(entry.path().string(), graph);
      graph_order.resize(graph.size());

      for (size_t id = 0; id < graph.size(); id++) {
        graph_order[graph[id]] = id;
      }
      std::string buf_file = fs::path(entry.path()).replace_filename("buf.bin").c_str();
      files_sys::bin::WriteSimple(buf_file, graph_order);

      char *loc_argv[4];
      loc_argv[0] = argv[0];
      loc_argv[1] = _strdup(file_vtk_out.c_str());
      loc_argv[2] = _strdup(buf_file.c_str());
      loc_argv[3] = _strdup(entry.path().filename().replace_extension().c_str());

      if (utils::AddScalarDataVtkFromBinFile(4, loc_argv)) {
        RETURN_ERR("Error when adding data\n");
      }

      fs::remove(buf_file);

      for (size_t i = 1; i < 4; i++) {
        free(loc_argv[i]);
      }
    }
  }

  return e_completion_success;
}
#endif // USE_VTK