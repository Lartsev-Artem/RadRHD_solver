#if !defined ILLUM_ADD_DIR_H
#define ILLUM_ADD_DIR_H

#include "solvers_struct.h"

namespace illum {

/**
 * @brief Пространство имён расчёта излучения по дополнительным направлениям для проецирования на картинную плоскость
 *
 */

namespace additional_direction {

/**
 * @brief Функция расчитывает переинтерполяцию сферы направлений на доп. направления и записывает результат в файл
 *
 * @param[in] base_address адрес с геометрией сферы направлений
 * @param[in] out_address  адрес для файлов переинтерполяции
 * @param[in] add_directions число дополнительных направлений
 * @param[in] grid_dir расчётная сфера направлений
 * @return int
 */
int MakeDirectionReInterpolation(const std::string &base_address, const std::string &out_address, const int add_directions, const grid_directions_t &grid_dir);

void SaveInterpolationScattering(const std::string &address_add_dir, const grid_directions_t &grid_dir, const grid_t &grid);
} // namespace additional_direction

} // namespace illum

#endif //! ILLUM_ADD_DIR_H