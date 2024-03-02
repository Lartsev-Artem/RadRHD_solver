#include "global_def.h"
#include "utils.h"

#include "geo_types.h"
#include "read_netgen.h"
#include "reader_txt.h"
#include "write_file_vtk.h"
#include "writer_txt.h"
#include <vector>

int FUNC_NAME(RenumNetgenByMetisToVtk)(int argc, char *argv[]) {
  std::string name_file_grid = "";
  std::string name_file_graph = "";
  std::string name_file_out = "";

  if (argc <= 3) {
    RETURN_ERR("Error input data!\n Input: path\\file_mesh_netgen, path\\file_graph, path\\output_file.vtk\n");
  } else {
    name_file_grid = argv[1];
    name_file_graph = argv[2];
    name_file_out = argv[3];
  }

  std::vector<int> graph;
  files_sys::txt::ReadData(name_file_graph, graph);

  std::vector<Vector3> point;
  std::vector<Eigen::Vector4i> cell;
  if (*(name_file_grid.end() - 1) == 'l') {
    utils::ReadNetgenMeshGrid(name_file_grid, point, cell);
  } else {
    utils::ReadNetgenGrid(name_file_grid, point, cell);
  }

  std::vector<std::pair<int, Eigen::Vector4i>> node_cells(cell.size());
  int i = 0;
  std::generate(node_cells.begin(), node_cells.end(),
                [&i, &graph, &cell] { i++; return std::make_pair(graph[i-1], cell[i-1]); });

  // std::cout << "def:\n";
  // for (size_t i = 0; i < node_cells.size(); i++) {
  //   std::cout << node_cells[i].first << ": " << node_cells[i].second[0] << node_cells[i].second[1] << node_cells[i].second[2] << node_cells[i].second[3] << "\n";
  // }

  std::sort(node_cells.begin(), node_cells.end(),
            [](std::pair<int, Eigen::Vector4i> l, std::pair<int, Eigen::Vector4i> r) { return l.first < r.first; });

  // std::cout << "\n\nres:\n";
  // for (size_t i = 0; i < node_cells.size(); i++) {
  //   std::cout << node_cells[i].first << ": " << node_cells[i].second[0] << node_cells[i].second[1] << node_cells[i].second[2] << node_cells[i].second[3] << "\n";
  // }

  i = 0;
  std::generate(cell.begin(), cell.end(),
                [&i, &node_cells] { i++; return node_cells[i-1].second; });

  i = 0;
  std::generate(graph.begin(), graph.end(),
                [&i, &node_cells] { i++; return node_cells[i-1].first; });

  files_sys::txt::WriteSimple(name_file_graph, graph);

  utils::WriteVtkFile(name_file_out, point, cell, 3);

  return e_completion_success;
}
