#ifdef MAKE_TRACE
#include "trace_main.h"
#ifdef TRANSFER_CELL_TO_FACE

#include "trace_nodes.h"
#include "intersections.h"
#include "trace_to_face.h"


int trace::RunTracesModule(TracerData& data, const Vector3& direction) {
    
  // make
  const std::vector<Face>& grid = data.grid;
  const std::vector<Matrix4> &vertexs = data.vertexs;
  const std::vector<Normals>& normals = data.normals;
  const std::vector<IntId>& neighbours = data.neighbours;
  
  const std::vector<BasePointTetra>& vec_x = data.vec_x;
  const std::vector<IntId>& sorted_graph = data.graph;
  const grid_t& geo_grid = data.geo_grid;

  // массивы записи в файл:
  std::vector<cell_local>& vec_x0 = data.vec_x0;    

  std::vector<IntId>& graph_bound_faces = data.graph_bound_faces;
  std::vector<graph_pair_t>& graph_cell_faces=data.graph_cell_faces;

  align_cell_local& X0 = data.X0;

  int num_cell;
  const int count_cells = vertexs.size();

  /*---------------------------------- далее FOR по направлениям----------------------------------*/
  for (int num_direction = 0; num_direction < 1; num_direction++)
  {      
    vec_x0.clear();
    vec_x0.reserve(CELL_SIZE * count_cells);

    graph_bound_faces.clear();
    graph_bound_faces.reserve(10000);
    graph_cell_faces.clear();
    graph_cell_faces.reserve(2 * geo_grid.cells.size());

    /*---------------------------------- далее FOR по ячейкам----------------------------------*/
    for (int h = 0; h < count_cells; ++h) {
      num_cell = sorted_graph[h];

      bits_flag_t face_state = 0;
      intersection::FindInAndOutFaces(direction, normals[num_cell], face_state);

      for (ShortId num_out_face = 0; num_out_face < CELL_SIZE; ++num_out_face) {
        if (CHECK_BIT(face_state, num_out_face) == e_face_type_out) // выходящие грани
        {
          GetLocNodes(num_cell, num_out_face, grid, vertexs[num_cell],
                      face_state, direction, normals[num_cell], neighbours,
                      vec_x[num_cell], vec_x0);

          graph_pair_t buf;
          buf.cell = num_cell;
          buf.loc_face = num_out_face;
          graph_cell_faces.push_back(buf); //это вместо graph
        }
        else if (neighbours[num_cell * CELL_SIZE + num_out_face] < 0) ///< сосед к текущей грани
        {
          graph_bound_faces.push_back(geo_grid.cells[num_cell].geo.id_faces[num_out_face]); //номера граней на сетке
        }
      }
    }
    /*---------------------------------- конец FOR по ячейкам----------------------------------*/
    
    X0.s.resize(vec_x0.size());
    for (size_t i = 0; i < vec_x0.size(); i++)    
      X0.s[i] = vec_x0[i].s;

    X0.in_face_id.resize(vec_x0.size() / NODE_SIZE);            
    face_loc_id_t val;
    for (size_t i = 0; i < vec_x0.size() / NODE_SIZE; i++) {
      val.a = vec_x0[i * NODE_SIZE + 0].in_face_id;
      val.b = vec_x0[i * NODE_SIZE + 1].in_face_id;
      val.c = vec_x0[i * NODE_SIZE + 2].in_face_id;      
      X0.in_face_id[i].bits = val.bits;      
    }      
  }
  /*---------------------------------- конец FOR по направлениям----------------------------------*/  
  return e_completion_success;
}
#endif  //! TRANSFER_CELL_TO_FACE
#endif  //! MAKE_TRACE