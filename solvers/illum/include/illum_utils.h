/**
 * @file illum_utils.h
 * @brief Файл включает функции общие для всех реализация алгоритма
 */

#if !defined ILLUM_UTILS_H && defined ILLUM && defined SOLVERS
#define ILLUM_UTILS_H

#include "global_types.h"
#include "solvers_struct.h"

/*! \addtogroup illum Модуль расчёта излучения
    \brief Модуль содержит функции для расчёта излучения в трех мерной области методом коротких характеристик
    @{
*/

/**
 * @brief Пространство имён модуля расчёта излучения
 *
 */
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

/**
 * @brief Функция возвращает среднее значение по ячейке
 *
 * @param[in] array указатель на начало данных длинны ::CELL_SIZE
 * @return среднее значение по ячейке
 */
inline Vector3 GetAverageByCell(const Vector3 *array) {
  Vector3 val = Vector3::Zero();
  for (int i = 0; i < CELL_SIZE; i++) {
    val += array[i];
  }
  return (val / CELL_SIZE);
}

/**
 * @brief Функция возвращает среднее значение по ячейке
 *
 * @param[in] array указатель на начало данных длинны ::CELL_SIZE
 * @return среднее значение по ячейке
 */
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

/**
 * @brief Функция возвращает значение на определяющей грани
 *
 * @param[in] num_in_face номер входящей грани(локальный)
 * @param[in] faces все грани
 * @param[in] cell текущая ячейка
 * @param[in] inter_coef коэффицинты интерполяции текущей грани
 * @return определяющее значение на входящей грани для текущего выходящего узла
 */
Type GetIllumeFromInFace(const int num_in_face, const int neigh_id, elem_t *cell, Vector3 &inter_coef);
} // namespace illum

#endif //! ILLUM_UTILS_H