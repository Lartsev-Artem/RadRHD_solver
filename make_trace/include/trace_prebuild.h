/**
 * @file trace_prebuild.h
 * @brief Функции предварительного расчёта геометрии сетки
 *
 */

#if !defined TRACE_PREBUILD_H
#define TRACE_PREBUILD_H

#include "json_struct.h"

namespace trace {

/**
 * @brief Функция строит матрицы перехода к локальному тетраэдру и
 * формирует список всех граней с их координатами вершин
 *
 * @param glb_files структура файлов (инициализированная)
 * @return int ::e_type_completion
 */
int PreBuild(const global_files_t &glb_files);

} // namespace trace
#endif //! TRACE_PREBUILD_H