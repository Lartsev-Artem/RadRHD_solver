/**
 * @file ray_tracing_build_plane.h
 * @brief Построение лучей, формирующих картинную плоскость
 *
 */
#ifndef RAY_TRACING_BUILD_PLANE_H
#define RAY_TRACING_BUILD_PLANE_H

#include "geo_types.h"

namespace ray_tracing {
/**
 * @brief Функция создает лучи для картинной плоскости на текущем кадре
 *
 * @param[in] num_frame номер кадра
 * @param[out] rays набор лучей на каждый пиксель
 */
void MakePlane(int num_frame, std::vector<Ray_t> &rays);
} // namespace ray_tracing
#endif //! RAY_TRACING_BUILD_PLANE_H