/**
 * @file radRHD_main.h
 * @brief Файл подключает модуль расчета радиационной газовой динамики в релятивистской постановке c
 * @version 0.1
 * @date 2023-10-31
 *
 */
#if !defined RADRHD_MAIN_H && defined RAD_RHD
#define RADRHD_MAIN_H

/*! \addtogroup rad_rhd Модуль расчета радиационной газовой динамики в релятивистской постановке
    @{
*/

/**
 * @brief Пространство имён модуля радиационной газовой динамики в релятивистской постановке
 *
 */
namespace rad_rhd {

int RadRHD_ConstRadStateTest();
} // namespace rad_rhd

#endif //! RADRHD_MAIN_H