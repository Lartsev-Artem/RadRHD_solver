/**
 * @file rhllc_init.h
 * @brief Инициализация hllc решателя
 *
 */
#if !defined RHLLC_INIT_H && defined SOLVERS
#define RHLLC_INIT_H

#include "solvers_struct.h"

/*! \addtogroup rhllc Модуль расчета газовой динамики в релятивистской постановке
    @{
*/

namespace rhllc {
/**
 * @brief Устанавливает настройки решателя по умолчанию.
 *
 * @param[out] hllc_set
 * @note доступа динамическая инициализация через json-config файл
 */
void SetCfgDefault(hllc_value_t &hllc_set);

/**
 * @brief Функция устанавливает начальное состояние газа
 *
 * @param file_init_value файл с предварительным газодинамическим состоянием
 * @param cells ячейки сетки
 * @return int  ::e_type_completion
 * @note в случае отсутствия файла инициализации будет установлено распределение по умолчанию
 */
int Init(std::string &file_init_value, std::vector<elem_t> &cells);
} // namespace rhllc
#endif //! RHLLC_INIT_H