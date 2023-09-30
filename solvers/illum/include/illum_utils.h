#if !defined ILLUM_UTILS_H && defined ILLUM && defined SOLVERS
#define ILLUM_UTILS_H

#include "global_types.h"
#include "solvers_struct.h"

namespace illum {

/**
 * @brief Функция возвращает среднее значение по ячейке
 *
 * @param[in] array указатель на начало данных длинны ::CELL_SIZE
 * @return среднее значение по ячейке
 */
inline double GetAverageByCell(const double *array) {
  double val = 0;
  for (int i = 0; i < CELL_SIZE; i++) {
    val += array[i];
  }
  return (val / CELL_SIZE);
}

inline Vector3 GetAverageByCell(const Vector3 *array) {
  Vector3 val = Vector3::Zero();
  for (int i = 0; i < CELL_SIZE; i++) {
    val += array[i];
  }
  return (val / CELL_SIZE);
}
inline Matrix3 GetAverageByCell(const Matrix3 *array) {
  Matrix3 val = Matrix3::Zero();
  for (int i = 0; i < CELL_SIZE; i++) {
    val += array[i];
  }
  return (val / CELL_SIZE);
}

/**
 * @brief Функция возвращает значение на границе
 *
 * @param[in] type_bound тип граничных условий  ::e_boundary_types_t
 * @return значение излучения
 */
Type BoundaryConditions(const int type_bound);

/**
 * @brief Функция возвращает текущее накопленное значение вдоль луча
 *
 * @param[in] x точка на сетке в которой рассчитывается излучение
 * @param[in] s расстояние пройденное лучом от начала излучения
 * @param[in] I_0 значение излучения на старте луча
 * @param[in] int_scattering интеграл рассеяния в ячейке
 * @param[in] cell ячейка с газодинамическими параметрами
 * @return значение излучения
 * @warning завязано на тип геометрии и конфигурацию решения ::e_grid_vtk_config_t
 */
Type GetIllum(const Vector3 x, const Type s, const Type I_0, const Type int_scattering, const elem_t &cell);

/**
 * @brief Функция переводит данные расчёта излучения в структуру решателя
 *
 * @note используется только при расчёте на cpu, для расчёта на видеокарте не нужна
 * @param num_dir номер направления
 * @param inter_coef коэффициенты интерполяции по данному направлению
 * @param grid сетка с излучением
 * @return возвращает норму ошибки на текущей итерации
 * @warning сейчас коэффициенты интерполяции == значения в узлах на гранях
 */
Type ReCalcIllum(const int num_dir, const std::vector<Vector3> &inter_coef, grid_t &grid);
} // namespace illum

#endif //! ILLUM_UTILS_H