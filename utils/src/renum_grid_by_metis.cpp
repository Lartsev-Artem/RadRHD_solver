#include "global_def.h"
#include "utils.h"

#include "geo_types.h"
#include "read_netgen.h"
#include "reader_txt.h"
#include "write_file_vtk.h"
#include "writer_txt.h"
#include <vector>

#ifdef USE_VTK
#include "reader_vtk.h"
#include "convert_vtk_geo.h"
/**
 * @brief Сортировка граничащих ячеек по краям подобластей
 * @warning Поведение если одна ячейка граничит сразу с тремя подобластями НЕОПРЕДЕЛЕННО!!!
 * @param name_file_vtk сетка
 * @param graph metis нумерация ячеек
 * @param cells ячейки сетки
 * @return int 
 */
static int sort_bound_cells(const std::string& name_file_vtk, const std::vector<int>& graph, std::vector<Eigen::Vector4i>& cells)
{ 
  vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  if(files_sys::vtk::Read(name_file_vtk,unstructured_grid))
  {
    printf("Error vtk reading\n");
    return e_completion_fail;
  }

  std::vector<int> neighbors;
  if (GetNeighborFace3D(unstructured_grid, neighbors))
  {
    printf("Error neighbor building\n");
    return e_completion_fail;
  }

  const int nodes = *std::max_element(graph.begin(), graph.end()) + 1;
  
  std::vector<int> size_nodes(nodes,0);
  for (size_t i = 0; i < nodes; i++)
  {
      size_nodes[i] = std::accumulate(graph.begin(), graph.end(),0 , [&i](int in,int a){return in+(int)(a==i);});      
  }

  int start = 0;
  int end = 0;
  for (int np = 0; np < nodes; np++)
  {
    start = end;
    end += size_nodes[np];
    
    auto neigh = neighbors.begin() + start*CELL_SIZE;

    auto head = cells.begin() + start;
    auto tail = cells.begin() + end-1;
    auto cur = head;
        
    for (; cur != tail; cur++, neigh+=CELL_SIZE)
    {
      for (int k = 0; k < CELL_SIZE; k++)
      {          
        if(neigh[k] >= 0) //не внешняя граница
        {
          int neigh_np = graph[neigh[k]/CELL_SIZE];             
          if(neigh_np < np) //соседний узел слева
          { 
            int h = CELL_SIZE*std::distance(cells.begin(), head);                
            std::swap_ranges(neighbors.begin()+h, neighbors.begin()+h + 4, neigh);
            std::swap(*head, *cur);
            ++head;
            break;
          }
          else if(neigh_np > np)//соседний узел справа
          { 
            int t = CELL_SIZE*std::distance(cells.begin(), tail);                
            std::swap_ranges(neighbors.begin()+t, neighbors.begin()+t +4, neigh);
            std::swap(*tail, *cur);
            --cur;
            --tail;
            neigh-=CELL_SIZE;
            break;
          }            
        }
      }      
    }         
  }

  return e_completion_success;
}
#endif //USE_VTK

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

#ifdef USE_VTK
  if(sort_bound_cells(name_file_out, graph, cell) == e_completion_success)
  {
    utils::WriteVtkFile(name_file_out, point, cell, 3); //пишем перенумерованную сетку
  }
#endif
  
  return e_completion_success;
}