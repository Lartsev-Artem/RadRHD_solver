/**
 * @file illum_utils.h
 * @brief Файл включает функции общие для всех реализация алгоритма
 */

#if !defined ILLUM_UTILS_H && defined ILLUM && defined SOLVERS
#define ILLUM_UTILS_H

#include "geo_types.h"
#include "global_types.h"
#include "ray_tracing_const.h"
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
 * @param[in] type_bound тип граничных условий  ::e_boundary_types_t (связан с геометрией сетки)
 * @param[in] type_obj   тип объекта, определяющего условие на границе ::e_ray_intersect_code
 * @param[in] inter_coef  коэффициенты интерполяции на потонциально определяющий границу грани
 * @return значение излучения
 */
Type BoundaryConditions(const IdType type_bound, const IntId type_obj = e_ray_intersect_none, const Vector3 &inter_coef = Vector3 ::Zero());

/**
 * @brief Функция возвращает текущее накопленное значение вдоль луча
 *
 * @param[in] x точка на сетке в которой рассчитывается излучение
 * @param[in] s расстояние пройденное лучом от начала излучения
 * @param[in] I_0 значение излучения на старте луча
 * @param[in] int_scattering интеграл рассеяния в ячейке
 * @param[inout] cell ячейка с газодинамическими параметрами
 * @return значение излучения
 * @warning завязано на тип геометрии и конфигурацию решения ::e_grid_vtk_config_t
 */
Type GetIllum(const Vector3 x, const Type s, const Type I_0, const Type int_scattering, elem_t &cell);

/**
 * @brief Функция переводит данные расчёта излучения в структуру решателя
 *
 * @note используется только при расчёте на cpu, для расчёта на видеокарте не нужна
 * @param[in] num_dir номер направления
 * @param[in] inter_coef коэффициенты интерполяции по данному направлению
 * @param[inout] grid сетка с излучением
 * @param[in] mpi_dir_shift - сдвиг в массиве излучения
 * @note при расчёте mpi не вся сетка расчитывается на узле, но вся собирается
 * @return возвращает норму ошибки на текущей итерации
 * @warning сейчас коэффициенты интерполяции == значения в узлах на гранях
 */
Type ReCalcIllum(const IdType num_dir, const std::vector<Vector3> &inter_coef, grid_t &grid, IdType mpi_dir_shift = 0);

#ifdef TRANSFER_CELL_TO_FACE

/**
 * @brief Функция переводит данные расчёта излучения в структуру решателя
 *
 * @param[in] num_dir номер направления
 * @param[in] inter_coef значения на гранях  по данному направлению
 * @param[inout] grid сетка с излучением
 * @param[in] mpi_dir_shift - сдвиг в массиве излучения
 * @return возвращает норму ошибки на текущей итерации
 */
Type ReCalcIllum(const IdType num_dir, const std::vector<Type> &inter_coef, grid_t &grid, IdType mpi_dir_shift = 0);
/**
 * @brief Функция возвращает значение на определяющей грани
 *
 * @param[in] faces все грани
 * @param[in] inter_coef коэффицинты интерполяции текущей грани
 * @return определяющее значение на входящей грани для текущего выходящего узла
 */
Type GetIllumeFromInFace(const IdType neigh_id, Vector3 &inter_coef
#ifdef INTERPOLATION_ON_FACES
                         ,
                         const Vector2 &x0 = Vector2::Zero()
#endif
);

namespace separate_gpu {
/**
 * @brief Функция переводит данные расчёта излучения в структуру решателя
 *
 * \note Пишет сразу в 2 массива. Глобальным с нумераций по ячейкам
 * и локальный массив на отправку по направлениям
 * @param[in] num_dir номер направления (локальный)
 * @param[in] inter_coef значения на гранях по данному направлению
 * @param[inout] grid сетка с излучением
 * @param[in] dir_disp - сдвиг по локальным направлениям
 * @return возвращает норму ошибки на текущей итерации
 */
Type ReCalcIllum(const IdType num_dir, const std::vector<Type> &inter_coef, grid_t &grid, const IdType dir_disp = 0);
} // namespace separate_gpu

/**
 * @brief
 *
 * @param[in] x точка на сетке в которой рассчитывается излучение(центр ячейки)
 * @param[in] int_scattering интеграл рассеяния в ячейке
 * @param[inout] cell ячейка с газодинамическими параметрами
 * @param[out] k коэффициент ослабления (alpha+betta)
 * @return правая часть : (alpha * Q + betta * S) / k;
 */
Type GetRhs(const Vector3 x, const Type int_scattering, elem_t &cell, Type &k);

/**
 * @brief Расчет среднего значения на грани
 *
 * @note  использует векторные регистры 256
 * @warning исключен предельный переход, который есть в скалярной версии
 * @param[in] I0 вектор значений излучения на старте луча
 * @param[in] s вектор расстояний пройденное лучом от начала излучения
 * @param[in] k коэффициент ослабления (alpha+betta)
 * @param[in] rhs правая часть : (alpha * Q + betta * S) / k;
 * @return Type излучение на грани
 */
Type GetIllum(const Type *I0, const Type *s, const Type k, const Type rhs);
#endif

} // namespace illum

#endif //! ILLUM_UTILS_H