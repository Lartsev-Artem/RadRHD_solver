/**
 * @file ray_tracing_main.h
 * @brief  Файл подключает модуль трассировки лучей и формирования картинной плоскости
 * @version 0.1
 * @date 2023-10-12
 *
 */
#ifndef RAY_TRACING_MAIN_H
#define RAY_TRACING_MAIN_H

/**
 * @brief Пространство имён трассировки лучей на картинную плоскость
 *
 */
namespace ray_tracing {

/**
 * @brief Функция запускает геометрическую сортировку
 * и формирует файлы матриц изображения
 *
 * @return int ::e_type_completion
 */
int RunRayTracing();
} // namespace ray_tracing
#endif //! RAY_TRACING_MAIN_H