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

int BuildSphereDirection(const std::string &sphere_vtk_grid, const std::string &out_address);

#endif

#endif //! MAKE_INTERNAL_FORMAT_H