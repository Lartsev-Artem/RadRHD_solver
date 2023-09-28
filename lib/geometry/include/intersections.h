/**
 * @file intersections.h
 * @brief Функции поиска пересечений геометрических объектов
 */
#ifndef INTERSECTIONS_H
#define INTERSECTIONS_H
#include "geo_types.h"
#include "global_types.h"

/**
 * @brief Пространство имён геометрических пересечений
 *
 */
namespace intersection {

/**
 * @brief Функция определяет точку пересечения луча с плоскостью, заданной 3 точками
 *
 * @param[in] face плоскость
 * @param[in] start_point начало луча
 * @param[in] direction направление
 * @param[out] result точка пересечения
 */
void IntersectionWithPlane(const Face &face, const Vector3 &start_point, const Vector3 &direction, Vector3 &result);

/**
 * @brief Функция проверяет пересечение треугольника и точки на плоскости
 *
 * @param[in] cell_face плоскость
 * @param[in] normal внешняя нормаль к плоскости грани
 * @param[in] XX проверочная точка
 * @return bool признак принадлежности
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
} // namespace intersection
#endif //! INTERSECTIONS_H