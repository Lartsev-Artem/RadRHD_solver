/**
 * @file ray_tracing_const.h
 * @brief Константы с описанием геометрии
 *
 */
#ifndef RAY_TRACING_CONST_H
#define RAY_TRACING_CONST_H

#include "geo_types.h"

/*! \addtogroup ray_tracing Модуль трассировки лучей, построения картинной плоскости и кривых блеска
    @{
*/

/**
 * @brief Код объекта пересечения
 * \note неотрицательные значения будут соответствовать ячейкам (или граням)
 *
 */
enum e_ray_intersect_code {
  e_ray_intersect_none = -1,   ///< нет пересечения
  e_ray_intersect_disk = -2,   ///< пересечение с диском
  e_ray_intersect_sphere = -3, ///< пересечение со сферой
  e_ray_intersect_rosh = -4    ///< пересечение с потенциалом роша
};

#ifdef __NVCC__
#define HOST_DEVICE __device__
#else
#define HOST_DEVICE
#endif

namespace ray_tracing {

HOST_DEVICE constexpr double PI = 3.1415926535897932384626433832795; ///<число пи

HOST_DEVICE constexpr Type k_center_x = 1, k_center_y = 0, k_center_z = 0; ///< координаты центр аккретора
const Vector3 k_center_sphere(k_center_x, k_center_y, k_center_z);         ///< центр аккретора
HOST_DEVICE constexpr Type k_radius_sphere = 0.01;                         ///< радиус аккретора
HOST_DEVICE constexpr Type k_internal_radius_disk = 0.01;                  ///< внутренний радиус около аккриционного диска
HOST_DEVICE constexpr Type k_external_radius_disk = 0.06;                  ///< внешний радиус около аккриционного диска

// параметры картинной плоскости
constexpr Type k_width_plane = 2.4;         ///< безразмерная ширина картинной плоскости
constexpr Type k_height_plane = 1.2;        ///< безразмерная высота картинной плоскости
constexpr int k_pixels_width = 200;         ///< число пикселей в ширину
constexpr int k_pixels_height = 100;        ///< число пикселей в высоту
constexpr Type k_height_above_center = 0.3; ///< высота над центром масс

const Vector3 k_center_of_mass(0.879518, 0, 0); ///< центр масс системы
constexpr int k_number_of_frame = 10;           ///< число позиций в окружности

} // namespace ray_tracing
#endif //! RAY_TRACING_CONST_H