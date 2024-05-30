/**
 * @file ray_tracing_main.h
 * @brief  Файл подключает модуль трассировки лучей и формирования картинной плоскости
 * @version 0.1
 * @date 2023-10-12
 *
 */
#ifndef RAY_TRACING_MAIN_H
#define RAY_TRACING_MAIN_H

/*! \addtogroup ray_tracing Модуль трассировки лучей, построения картинной плоскости и кривых блеска
    @{
*/

#include <string>
/**
 * @brief Пространство имён трассировки лучей на картинную плоскость
 *
 */
namespace ray_tracing {

/**
 * @brief Функция запускает геометрическую сортировку
 * и формирует файлы матриц изображения
 *
 * @param[in] file_energy файл с рассчитанной на сетке энергией
 * @return int ::e_type_completion
 */
int RunRayTracing(const std::string &file_energy);

/**
 * @brief Функция находит номера ячеек которые участвуют в проецировании на картинную плоскость
 * @note Параметры проекции задаются вручную до компиляции
 * @return int ::e_type_completion
 */
int FindIntersections();

/**
 * @brief Функция находит номера направлений, отраженные от внутренних стенок
 * @return int ::e_type_completion
 */
int FindReflectionDirections();
} // namespace ray_tracing
#endif //! RAY_TRACING_MAIN_H