/**
 * @file cuda_ray_geo_intersect.h
 * @brief Функции поиска пересечений с геометрическими объектами
 *
 */
#if !defined CUDA_RAY_GEO_INTERSECT_H && defined USE_CUDA
#define CUDA_RAY_GEO_INTERSECT_H

#include "cuda_ray_struct.h"

/*! \addtogroup ray_tracing Модуль трассировки лучей, построения картинной плоскости и кривых блеска
    @{
*/

namespace cuda::ray_tracing {

/**
 * @brief Расчёт параметра пересечния луча с поверхностью потенциала роша
 *
 * @param[in] ray луч
 * @param[in] t параметр вдоль луча?
 * @return значение ...
 */
__device__ Type Rosh(const Ray &ray, Type t);

/**
 * @brief Функция ищет корни уравнения на потенциала роша методом бисекции
 *
 * @param[in] n число разбиений
 * @param[in] a начало отрезка
 * @param[in] b конец отрезка
 * @param[in] ray луч
 * @param[out] roots корни, отсортированные по возрастанию. max_size=6
 * @return число корней
 * @warning массив для корней должен быть не менее 6
 */
__device__ int FindRoshRoots(const int n, const Type a, const Type b, const Ray &ray, Type *roots);

/**
 * @brief Возвращает признак пересечения луча с поверхностью потенциала роша
 *
 * @param[in] ray луч
 * @param[out] dist расстояние до пересечения
 * @return ::e_ray_intersect_code
 */
__device__ int GetIntersectionWithRosh(const Ray &ray, Type *dist);

/**
 * @brief Возвращает признак пересечения с диском или аккретором
 *
 * @param[in] ray луч
 * @return ::e_ray_intersect_code
 */
__device__ int GetIntersectionWithSphereODisk(const Ray &ray);

/**
 * @brief Пересечение луча с плоскость заданной 3 точками
 *
 * @param[in] face точки плоскости
 * @param[in] ray луч
 * @param[out] result точка пересечения
 */
__device__ void IntersectionWithPlane(const Face &face, const Ray &ray, Vector3 &result);

} // namespace cuda::ray_tracing
#endif //! CUDA_RAY_GEO_INTERSECT_H