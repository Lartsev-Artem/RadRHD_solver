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
 * @brief Возвращает отражённое направление
 *
 * @param[in] dir исходное направление
 * @param[in] n нормаль к поверхности (нормированная!)
 * @return направление отражения
 */
static inline Vector3 reflection_dir(const Vector3 &dir, const Vector3 &n)
{
  return (dir - 2.0 * n * (dir.dot(n))).normalized();
}

Vector3 Spherical2Cartesian(const Vector3 &rpt);
Vector3 Cartesian2Spherical(const Vector3 &xyz);

#endif //! LINEAR_ALG_H
