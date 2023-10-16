/**
 * @file hllc_main.h
 * @brief Файл подключает модуль расчета газовой динамики
 * @version 0.1
 * @date 2023-10-16
 *
 */
#if !defined HLLC_MAIN_H && defined SOLVERS && defined HLLC
#define HLLC_MAIN_H

/*! \addtogroup hllc Модуль расчета газовой динамики
    @{
*/

/**
 * @brief Пространство имён модуля газовой динамики
 *
 */
namespace hllc {

int RunHllcModule();

} // namespace hllc

#endif //! HLLC_MAIN_H