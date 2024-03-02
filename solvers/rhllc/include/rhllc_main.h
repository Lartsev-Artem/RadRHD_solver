/**
 * @file rhllc_main.h
 * @brief Файл подключает модуль расчета газовой динамики в релятивистской постановке
 * @version 0.1
 * @date 2023-10-07
 *
 */
#if !defined RHLLC_MAIN_H && defined SOLVERS
#define RHLLC_MAIN_H

/*! \addtogroup rhllc Модуль расчета газовой динамики в релятивистской постановке
    @{
*/

/**
 * @brief Пространство имён модуля газовой динамики в релятивистской постановке
 *
 */
namespace rhllc {
int RunRhllcModule();
int RunRhllcMpiModule();
} // namespace rhllc

#endif //! RHLLC_MAIN_H