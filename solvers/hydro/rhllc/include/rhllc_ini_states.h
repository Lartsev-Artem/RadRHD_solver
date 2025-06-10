/**
 * @file rhllc_ini_states.h
 * @author Artem
 * @brief Функции инициализации начального состояния для газовой динамики
 * @version 0.1
 * @date 2025-06-08
 *
 * @copyright Copyright (c) 2025
 *
 */
#if !defined RHLLC_INI_STATES_H && defined RHLLC
#define RHLLC_INI_STATES_H

#include "solvers_struct.h"

/*! \addtogroup rhllc Модуль расчета газовой динамики в релятивистской постановке
    @{
*/
namespace rrhd
{
    namespace rhllc
    {
        namespace ini
        {
            flux_t Soda(const Vector3 &x);
            flux_t Jet(const Vector3 &x);
            flux_t Uniform(const Vector3 &x);
        }; // namespace ini

        /**
         * @brief Функция устанавливает начальное состояние газа
         *
         * @param file_init_value файл с предварительным газодинамическим состоянием
         * @param ini_func функция инициализации
         * @param grid  сетка
         * @return int  ::e_type_completion
         * @note в случае отсутствия файла инициализации будет установлено распределение по умолчанию
         */
        int Init(const std::string &file_init_value,
                 const std::function<flux_t(const Vector3 &)> ini_func,
                 grid_t &grid);
    };
}; // namespace rrhd
#endif // RHLLC_INI_STATES_H