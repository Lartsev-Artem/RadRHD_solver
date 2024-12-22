/**
 * @file ray_tracing_build_plane.h
 * @brief Построение лучей, формирующих картинную плоскость
 *
 */
#ifndef RAY_TRACING_BUILD_PLANE_H
#define RAY_TRACING_BUILD_PLANE_H

#include "geo_types.h"
#include "ray_tracing_const.h"

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

/**
 * @brief Функция создает лучи для картинной плоскости на основе луча, выпущенного из центра плоскости
 *
 * @param[in] center_ray центральный луч
 * @param[out] rays набор лучей на каждый пиксель
 */
void MakeRays(const Ray_t &center_ray, std::vector<Ray_t> &rays);


/**
 * @brief Функция создает лучи для картинной плоскости выпущенные из положения наблюдателя и проходящие через пиксель
 *
 * @param[in] plane_params параметры картинной плоскости
 * @param[in] plane_orig центр картинной плоскости
 * @param[in] observer положение наблюдателя
 * @param[out] rays набор лучей на каждый пиксель
 */
void MakeRays(const PlaneParams& plane_params, const Vector3 &plane_orig, const Vector3& observer,std::vector<Ray_t> &rays);

#ifdef USE_VTK
/**
 * @brief Функция создает пустую сетку для картинной плоскости
 *
 * @param[in] x_cnt число пикселей в ширину
 * @param[in] y_cnt число пикселей в высоту
 * @param[out] image_plane сетка
 * @param[in] width_plane ширина плоскости(безразмерная)
 * @param[in] height_plane высота плоскости(безразмерная)
 */
void MakeVtkPlane(const int x_cnt, const int y_cnt, vtkSmartPointer<vtkUnstructuredGrid> &image_plane,
                  const Type width_plane =k_width_plane, const Type height_plane=k_height_plane);

void MakeVtkPlane(const PlaneParams& plane, vtkSmartPointer<vtkUnstructuredGrid> &image_plane);
/**
 * @brief Функция перестраивает бинарные файлы в формат визуализации vtk (получает серию картинных плоскостей)
 *
 * @param[in] plane_cfg параметры плоскости  
 * @param[in] number_of_planes число плоскостей
 * @param[in] files_plane основное название серии бинарных файлов с картинной плоскостью (в формате path/file[i.bin])
 * @return int ::e_type_completion
 */
int BuildVtkFromBin(const PlaneParams& plane_cfg, const int number_of_planes, const std::string &files_plane);
#endif
} // namespace ray_tracing
#endif //! RAY_TRACING_BUILD_PLANE_H