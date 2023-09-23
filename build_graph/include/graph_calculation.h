#if !defined GRAPH_CALCULATION_H && BUILD_GRAPH
#define GRAPH_CALCULATION_H

#include "geo_types.h"
#include <set>

// #include <bitset>
// #include <list>
// #include <map>
// #include <set>
// #include <vector>

// #include "geo_types.h"
// #include "global_types.h"
// #include "graph_config.h"

// #ifdef USE_OMP
// #include <omp.h>
// #endif

namespace graph {

//=======================================OMP=======================================
#ifdef USE_OMP

#ifdef GRID_WITH_INNER_BOUNDARY
int FindCurCellWithHole(const std::vector<IntId> &next_step_el, const std::vector<IntId> &count_in_face, std::vector<IntId> &count_knew_face,
                        std::list<IntId> &cur_el, const std::set<IntId> &inner_part, std::set<IntId> &outter_part,
                        const std::map<IntId, FaceCell> &inner_cells, const std::vector<IntId> &all_pairs_id, const Vector3 &direction, const std::vector<Normals> &normals);
#else
int FindCurCell(const std::vector<IntId> &next_step_el, const std::vector<IntId> &count_in_face, const std::vector<IntId> &count_knew_face,
                std::list<IntId> &cur_el);
#endif

int FindNumberOfAllInnerFaceAndKnew(const Vector3 &dir, const std::vector<Normals> &normals, const std::vector<State> &faces_state,
                                    std::vector<IntId> &count_in_face, std::vector<IntId> &count_knew_face, std::vector<IntId> &next_step_el);

int NewStep(const std::vector<IntId> &all_pairs_id, const std::vector<IntId> &count_in_face, std::vector<IntId> &count_knew_face, std::list<IntId> &cur_el,
            std::vector<IntId> &next_step_el);

#else

#ifdef GRID_WITH_INNER_BOUNDARY
int FindCurCellWithHole(const std::set<IntId> &next_step_el, const std::vector<IntId> &count_in_face, std::vector<IntId> &count_knew_face,
                        std::vector<IntId> &cur_el,
                        const std::set<IntId> &inner_part, std::set<IntId> &outter_part, const std::map<IntId, FaceCell> &inner_cells,
                        const std::vector<IntId> &all_pairs_id, const Vector3 &direction, const std::vector<Normals> &normals);
#else
int FindCurCell(const std::set<IntId> &next_step_el, const std::vector<IntId> &count_in_face, const std::vector<IntId> &count_knew_face, std::vector<IntId> &cur_el);
#endif

int FindNumberOfAllInnerFaceAndKnew(const Vector3 &dir, const std::vector<Normals> &normals, const std::vector<State> &faces_state,
                                    std::vector<IntId> &count_in_face, std::vector<IntId> &count_knew_face, std::set<IntId> &next_step_el);

int NewStep(const std::vector<IntId> &all_pairs_id, const std::vector<IntId> &count_in_face, std::vector<IntId> &count_knew_face, const std::vector<IntId> &cur_el,
            std::set<IntId> &next_step_el);

#endif // USE_OMP

} // namespace graph
#endif //! GRAPH_CALCULATION_H
