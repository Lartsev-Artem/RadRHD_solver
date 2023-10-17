/**
 * @file linear_alg.h
 * @brief Функции линейной алгебры
 *
 */
#ifndef LINEAR_ALG_H
#define LINEAR_ALG_H

#include "geo_types.h"

/**
 * @brief Расчёт матрицы поворота к ячейке (для расчёта потоков по нормали)
 *
 * @note В полном виде матрица 5x5 с элементами a_11=a_55=1.
 * Т.к. давление и плотность скалярные величины. (т.е. поворачиваем только скорость)
 * @param[in] n нормаль
 * @param[out] T матрица поворота
 */
void GetRotationMatrix(const Vector3 &n, Matrix3 &T);

/**
 * @brief Возвращает центр треугольника
 *
 * @param[in] a вершина треугольника
 * @param[in] b вершина треугольника
 * @param[in] c вершина треугольника
 * @return точка - центр
 */
inline Vector3 GetCenterTriangle(const Vector3 &a, const Vector3 &b, const Vector3 &c) {
  return (a + b + c) / .3;
}
/**
 * @brief  Возвращает центр треугольника
 *
 * @param f - треугольник, заданный вершинами
 * @return точка - центр
 */
inline Vector3 GetCenterTriangle(const Face f) {
  return (f.A + f.B + f.C) / .3;
}

#endif //! LINEAR_ALG_H