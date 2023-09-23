#if !defined GRAPH_READ_WRITE && defined BUILD_GRAPH
#define GRAPH_READ_WRITE

#include "graph_config.h"

#include <map>
#include <memory>
#include <string>

#include "geo_types.h"
#include "global_types.h"

namespace graph {

#ifdef USE_STRANGE_FUNCTION
int ReadInitBoundarySet(const std::string name_file_boundary, std::set<IntId> &boundary_cells);
int ReadInnerBoundary(const std::string name_file_boundary_inner, std::set<IntId> &id_inner_boundary_face);
int ReadInnerCellOfSphere(const std::string name_file_inner_sphere, std::vector<Face> &inner_faces);
int ReadInnerCellOfSphereAndId(const std::string name_file_face_and_id, std::map<IntId, Face> &inner_faces);
#endif

} // namespace graph
#endif // BUILD_GRAPH_READ_WRITE