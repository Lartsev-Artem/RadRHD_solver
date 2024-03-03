#if !defined SPECTRUM_UTILS_H && defined SPECTRUM

#include "solvers_struct.h"

namespace illum {
namespace spectrum {

/**
 * @brief Функция возвращает полное излучение по указанному направлению нормированное на расстояние
 *
 * @param[in] num_dir номер направления
 * @param[in] grid сетка
 * @return излучение по всей сетки
 */
Type get_full_illum(const IdType num_dir, const grid_t &grid);

/**
 * @brief Инициализирует состояние на основе данных с диска
 *
 * @param[in] num индекс решения
 * @param[inout] grid сетка
 * @return код ошибки
 */
int InitPhysState(const int num, grid_t &grid);
}; // namespace spectrum
}; // namespace illum

#endif //! SPECTRUM_UTILS_H