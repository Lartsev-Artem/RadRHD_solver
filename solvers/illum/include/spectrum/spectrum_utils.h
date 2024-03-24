#if !defined SPECTRUM_UTILS_H && defined SPECTRUM
#define SPECTRUM_UTILS_H

#include "solvers_struct.h"

#ifdef DEBUG
//#define LOG_SPECTRUM
#endif

#ifdef LOG_SPECTRUM
extern bool log_enable;
#define log_spectrum(...) WRITE_LOG(__VA_ARGS__)
#else
#define log_spectrum(...)
#endif

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

/**
 * @brief
 *
 * @param[in] x точка на сетке в которой рассчитывается излучение(центр ячейки)
 * @param[in] int_scattering интеграл рассеяния в ячейке
 * @param[inout] cell ячейка с газодинамическими параметрами
 * @param[out] k коэффициент ослабления (alpha+betta)
 * @param[in] frq0 левая граница частоты
 * @param[in] frq1 правая граница частоты
 * @return правая часть : (alpha * Q +  S) / k;
 */
Type GetRhsOpt(const Vector3 x, const Type S, elem_t &cell, Type &k, Type frq0, Type frq1);
}; // namespace illum

#endif //! SPECTRUM_UTILS_H