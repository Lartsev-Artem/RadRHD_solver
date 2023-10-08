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

#endif //! LINEAR_ALG_H