/**
 * @file graph_static_value.cpp
 * @brief Файл содержит глобальные переменные модуля вычислителя графов
 *
 */
#ifdef BUILD_GRAPH
#include "graph_struct.h"

namespace graph {

std::vector<boundary_trace_t> bound_trace;

uint64_t id_try_size;
uint64_t dist_try_size;
uint64_t x_try_size;

} // namespace graph
#endif