/**
 * @file ray_tracing_build_plane.h
 * @brief Построение лучей, формирующих картинную плоскость
 *
 */
#ifndef RAY_TRACING_BUILD_PLANE_H
#define RAY_TRACING_BUILD_PLANE_H

#include "geo_types.h"

#ifdef USE_VTK
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#endif

/*! \addtogroup ray_tracing Модуль трассировки лучей, построения картинной плоскости и кривых блеска
    @{
*/

namespace ray_tracing {
/**
 * @brief Функция создает лучи для картинной плоскости на текущем кадре
 *
 * @param[in] num_frame номер кадра
 * @param[out] rays набор лучей на каждый пиксель
 */
void MakeRays(int num_frame, std::vector<Ray_t> &rays);

#ifdef USE_VTK
/**
 * @brief Функция создает пустую сетку для картинной плоскости
 *
 * @param[out] image_plane сетка
 */
void MakeVtkPlane(vtkSmartPointer<vtkUnstructuredGrid> &image_plane);

/**
 * @brief Функция перестраивает бинарные файлы в формат визуализации vtk (получает серию картинных плоскостей)
 *
 * @param[in] number_of_planes число плоскостей
 * @param[in] files_plane основное название серии бинарных файлов с картинной плоскостью (в формате path/file[i.bin])
 * @return int ::e_type_completion
 */
int BuildVtkFromBin(const int number_of_planes, const std::string &files_plane);
#endif
} // namespace ray_tracing
#endif //! RAY_TRACING_BUILD_PLANE_H