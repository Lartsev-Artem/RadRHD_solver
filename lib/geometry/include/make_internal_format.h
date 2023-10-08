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

#endif

#endif //! MAKE_INTERNAL_FORMAT_H