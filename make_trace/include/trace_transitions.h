#if !defined TRACE_CALC_H && defined MAKE_TRACE
#define TRACE_CALC_H

#include "geo_types.h"
#include "global_types.h"

namespace trace {

/**
 * @brief Функция преобразует точку в глобальных координатах в локальные координаты тетраэдра
 *
 * @param[in] vertex_tetra матрица перехода [X;Y;Z;1]
 * @param[in] global_coord точка в глобальных координатах
 * @param[out] local_coord точка в локальных координатах
 */
void FromGlobalToLocalTetra(const Eigen::Matrix4d &vertex_tetra, const Vector3 &global_coord, Vector3 &local_coord);

/**
 * @brief Функция преобразует точку в локальных координатах в глобальные координаты тетраэдра
 *
 * @param[in] vertex_tetra матрица перехода [X;Y;Z;1]
 * @param[in] local_coord точка в локальных координатах
 * @param[out] global_coord точка в глобальных координатах
 */
void FromLocalToGlobalTetra(const Eigen::Matrix4d &vertex_tetra, const Vector3 &local_coord, Vector3 &global_coord);

/**
 * @brief Функция преобразует точку из координат тетраэдра в координаты плоскости-грани тетраэдра
   @details по формуле:  \f$X_{plane} = T * (X_{tetra} - O)\f$, где \f$T\f$ --- матрица перехода, O - начало системы координат
 * @param[in] transform_matrix матрица перехода
 * @param[in] start_point точка начала координат
 * @param[in] tetra_coord точка в координатах тетраэдра
 * @param[out] plane_coord точка в координатах плоскости
 */
inline void FromTetraToPlane(const Matrix3 &transform_matrix, const Vector3 &start_point, const Vector3 &tetra_coord, Vector3 &plane_coord);

/**
 * @brief Функция преобразует точку из координат плоскости-грани в координаты тетраэдра
 *
 * @details по формуле:  \f$X_{tetra} = T^{-1} * X_{tetra} + O\f$, где \f$T\f$ --- матрица перехода  O - начало системы координат
 * tetra_coord = inverse_transform_matrix * plane_coord + start_point;
 * @param[in] inverse_transform_matrix обратная матрица перехода
 * @param[in] start_point точка начала координат
 * @param[in] plane_coord точка в координатах плоскости
 * @param[out] tetra_coord точка в координатах тетраэдра
 */
inline void FromPlaneToTetra(const Matrix3 &inverse_transform_matrix, const Vector3 &start_point, const Vector3 &plane_coord, Vector3 &tetra_coord);

/**
 * @brief Функция возвращает коэффициенты интерполяции для представления решения в виде   ax+by+c
 *
 * @details по формуле:  по формуле: \f$y=A^{-1}*x\f$,  где А - матрица узлов, x- вектор значений
 * @param[in] interpolation_nodes постояные узлы интерполяции (формат (в координатах стандартного тетраэдра) {x_i, y_i, 1})
 * @param[out] function_value значение функции в трех точках
 * @return Vector3 коэффициенты
 * \warning эти функции не используются. Описание требует проверки!!!
 */
inline Vector3 GetInterpolationCoef(const Matrix3 &interpolation_nodes, const Vector3 &function_value);

/**
 * @brief Функция возвращает значения функции по коэффициентам интерполяции ax+by+c
 *
 * @details по формуле: \f$y=A*x\f$,  где А - матрица узлов, x- вектор коэффициентов
 * @param[in] interpolation_nodes постояные узлы интерполяции (формат (в координатах стандартного тетраэдра) {x_i, y_i, 1})
 * @param[out] function_value значение функции в трех точках
 * @return Vector3 коэффициенты
 * \warning эти функции не используются. Описание требует проверки!!!
 */
inline Vector3 GetInterpolationCoefInverse(const Matrix3 &interpolation_nodes, const Vector3 &function_value);

} // namespace trace
#endif //! TRACE_CALC_H