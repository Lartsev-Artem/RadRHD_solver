/**
 * @file make_internal_format.h
 * @brief Функция разбиения vtk сетки на бинарные файлы
 */
#ifndef MAKE_INTERNAL_FORMAT_H
#define MAKE_INTERNAL_FORMAT_H

#include "prj_config.h"
#ifdef USE_VTK

#include "json_struct.h"

int BuildDataFromVTK(const global_files_t &glb_files);

/**
 * @brief Перевод vtk поверхности в формат сферы направлений
 *
 * @param[in] sphere_vtk_grid - сетка
 * @param[in] out_address - адрес для сборки
 * @param[in] extended - расширенный режим геометрии, для условий отражения
 * @return int
 */
int BuildSphereDirection(const std::string &sphere_vtk_grid, const std::string &out_address, bool extended = false);

#endif

#endif //! MAKE_INTERNAL_FORMAT_H