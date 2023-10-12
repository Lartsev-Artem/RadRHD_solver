/**
 * @file cuda_ray_interface.h
 * @brief Файл хранит методы, через которые можно обращаться из методом
 * для gcc компилятора
 *
 */
#if !defined CUDA_RAY_INTERFACE_H && defined USE_CUDA
#define CUDA_RAY_INTERFACE_H

#include "geo_types.h"

namespace cuda {

/**
 * @brief Пространство имен функций трассировки лучей исполняемых на видеокарте
 *
 */
namespace ray_tracing {

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

/**
 * @brief Пространство имён интерфейсных функция модуля cuda::ray_tracing
 *
 */
namespace interface {

/**
 * @brief Поиск пересечений с заданными лучами.
 * @details Функция копирует rays_host на карту и intersections с карты.
 * Определяется ближайшее пересечение от источника
 * @param[in] rays_host лучи формирующие картинную плоскость
 * @param[out] intersections код пересечения ::e_ray_intersect_code (с номером грани)
 * @return int ::e_type_completion
 * @warning метод не оптимизирован для поиска по всей сетки. Применяется для пересечения с границей
 */
int StartTracing(const std::vector<Ray_t> &rays_host, std::vector<IntId> &intersections);

/**
 * @brief Инициализация видеокарты
 *
 * @param[in] faces_host грани сетки с которыми проверяем пересечения(с копированием на карту)
 * @param[in] rays_host структура матрицы лучей для картинной плоскости(только выделение памяти)
 */
void InitDevice(const std::vector<FaceCell> &faces_host, const std::vector<Ray_t> &rays_host);

/**
 * @brief удаление структур на видеокарте
 *
 */
void ClearDevice();

} // namespace interface
} // namespace ray_tracing
} // namespace cuda
#endif //! CUDA_RAY_INTERFACE_H