/**
 * @file graph_static_value.cpp
 * @brief Файл содержит глобальные переменные модуля вычислителя графов
 *
 */
#ifdef BUILD_GRAPH
#include "graph_struct.h"

namespace graph {

std::vector<int> id_try_surface;
std::vector<double> dist_try_surface;
std::vector<Eigen::Vector3d> x_try_surface;
try_solve_t buf_try;

uint64_t id_try_size;
uint64_t dist_try_size;
uint64_t x_try_size;

} // namespace graph
#endif