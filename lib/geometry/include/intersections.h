/**
 * @file intersections.h
 * @brief Функции поиска пересечений геометрических объектов
 */
#ifndef INTERSECTIONS_H
#define INTERSECTIONS_H
#include "geo_types.h"
#include "global_types.h"

#ifdef USE_VTK
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#endif

/**
 * @brief Пространство имён геометрических пересечений
 *
 */
namespace intersection {

/**
 * @brief Построение локального базиса плоскости
 *
 * @details по начальной точке и нормали строит локальный базис картинной плоскости (vec1, vec2).
    нормаль дана.
    задаем второй вектор ортогонально, по св-ву скалярного умножения
    третий вектор из векторного произведения
 *
 * @param[in] normal нормаль к плоскости
 * @param[out] basis новый базис (тройка векторов)
 */
void SetBasis(const Vector3 &normal, Matrix3 &basis);

/**
 * @brief Функция определяет точку пересечения луча с плоскостью, заданной 3 точками
 *
 * @param[in] face плоскость
 * @param[in] start_point начало луча
 * @param[in] direction направление
 * @param[out] result точка пересечения
 * @warning Устаревшие функции. Используйте ::RayIntersectsTriangle
 */
void IntersectionWithPlane(const Face &face, const Vector3 &start_point, const Vector3 &direction, Vector3 &result);

/**
 * @brief Функция проверяет пересечение треугольника и точки на плоскости
 *
 * @param[in] cell_face плоскость
 * @param[in] normal внешняя нормаль к плоскости грани
 * @param[in] XX проверочная точка
 * @return bool признак принадлежности
 * @warning Устаревшие функции. Используйте ::RayIntersectsTriangle
 */
bool InTriangle(const Face &cell_face, const Vector3 &normal, const Vector3 &XX);

/**
 * @brief Определяет типы граней ячейки по отношению к направлению
 *
 * @param direction[in] направление
 * @param normals_cell[in] нормали к ячейке
 * @param face_type[out] тип ячейки ::e_face_in_out_type_t
 */
void FindInAndOutFaces(const Vector3 &direction, const Normals &normals_cell, bits_flag_t &face_type);

/**
 * @brief Поиск пересечений луча с диско сферой вне расчетной области
 *
 * @param X0
 * @param n
 * @param res
 * @warning не использовать. Нужна проверка
 */
int IntersectionWithPlaneDisk(const Vector3 &X0, const Vector3 &n, Vector3 &res);

/**
 * @brief Функция поиска пересечения луча с треугольником
 *
 * @details Алгоритм Моллера — Трумбора для работы которого не требуется предварительное вычисление уравнения плоскости, содержащей треугольник.
 * @param[in] ray_orig начало луча
 * @param[in] ray_dir направление луча
 * @param[in] triangle треугольник
 * @param[out] intersection точка пересечения
 * @return Type расстояние до пересечения или -1
 */
Type RayIntersectsTriangle(const Vector3 &ray_orig, const Vector3 &ray_dir, const Face &triangle, Vector3 &intersection);

#ifdef USE_VTK
/**
 * @brief Функция возвращает неупорядоченные номера ячеек, пересекающихся лучом
 *
 * @param[in] ray  луч
 * @param[in] grid сетка
 * @param[out] id_cells id ячеек
 */
void GetIntersectionCellId(const Ray_t &ray, const vtkSmartPointer<vtkUnstructuredGrid> &grid, std::vector<IntId> &id_cells);
#endif

} // namespace intersection
#endif //! INTERSECTIONS_H