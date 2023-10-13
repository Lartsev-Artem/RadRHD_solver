/**
 * @file ray_tracing_calc_illum.h
 * @brief Расчёт энергии на картинной плоскости, формирование файла картинной плоскости, подсчёт кривой блеска
 *
 */
#ifndef RAY_TRACING_CALC_ILLUM_H
#define RAY_TRACING_CALC_ILLUM_H

#include <string>

/*! \addtogroup ray_tracing Модуль трассировки лучей, построения картинной плоскости и кривых блеска
    @{
*/

namespace ray_tracing {

#ifdef USE_VTK
/**
 * @brief Функция заполняет шаблон картиной плоскости и считает кривую блеска
 * @note требуется поддержка VTK
 *
 * @param[in] file_energy файл с распределением плотности энергии излучения на исходной сетки
 * @return int ::e_type_completion
 * @warning предварительно должна быть вызвана трассировка лучей и заполнены файловые структуры
 */
int MakeEnergyAndCurve(const std::string &file_energy);
#endif
} // namespace ray_tracing

#endif //! RAY_TRACING_CALC_ILLUM_H