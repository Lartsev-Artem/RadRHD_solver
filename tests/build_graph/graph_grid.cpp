/**
 * @file graph_grid.cpp
 * @brief Тест для модуля построителя графов
 * @version 0.1
 * @date 2023-09-26
 *
 * \details тест из файла vtk и сферы направлений строит граф обхода и создает новую сетку
 * со скалярными полями где значение --- позиция ячейки в графе
 *
 * \note требуется поддержка VTK и собранный build утилит
 *
 */

#include "build_internal_format.h"
#include "graph_main.h"

#include "global_value.h"
#include "reader_bin.h"
#include "writer_bin.h"
#ifdef USE_VTK
int RunTestGraph() {

  glb_files.base_address = "./";
  glb_files.graph_address = "./graph/graph";
  glb_files.name_file_sphere_direction = ".txt";
  glb_files.name_file_vtk = ".vtk";
  glb_files.Build();

  BuildDataFromVTK(glb_files);

  for (size_t i = 0; i < 126; i++) {
    files_sys::bin::ReadSimple(glb_files.graph_address + std::to_string(i) + ".bin")
  }

  // for(all files)
  set_bin_scalar_to_vtk("vtk", "graph.bin", "vtk_out", "graph_i");
  // SetScalarDataVtkFromFile(argc, argv);

  graph test
      DataArray->SetNumberOfTuples(n);
  for (size_t i = 0; i < n; i++) {
    DataArray->SetTuple1(vector_data[i], i);
  }

  graph::RunGraphModule();

  return 0;
}
#endif // USE_VTK