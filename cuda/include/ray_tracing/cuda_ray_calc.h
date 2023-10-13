/**
 * @file cuda_ray_calc.h
 * @brief Методы трассировки
 *
 */
#if !defined CUDA_RAY_CALC_H && defined USE_CUDA
#define CUDA_RAY_CALC_H

#include "cuda_ray_struct.h"

/*! \addtogroup ray_tracing Модуль трассировки лучей, построения картинной плоскости и кривых блеска
    @{
*/

namespace cuda::ray_tracing {

/**
 * @brief Функция поиска пересечения луча с треугольником
 *
 * @details Алгоритм Моллера — Трумбора для работы которого не требуется предварительное вычисление уравнения плоскости, содержащей треугольник.
 * @param[in] ray  луч
 * @param[in] triangle сетка с треугольниками
 * @param[out] intersection точка пересечения
 * @return расстояние до пересечения. (при отсутствии -1)
 */
__device__ Type RayIntersectsTriangle(const Ray &ray, const Face &triangle, Vector3 &intersection);

/**
 * @brief Поиск всех пересечений с матрицей лучей
 *
 * @param[in] M число лучей
 * @param[in] rays лучи
 * @param[in] N число треугольников
 * @param[in] triangles треугольники
 * @param[out] intersections коды пересечений для каждого луча
 * @return __global__ use cudaGetLastError
 */
__global__ void RayTracing(const int M, const Ray *rays, const int N, const Face *triangles, int *intersections);

} // namespace cuda::ray_tracing
#endif //! CUDA_RAY_CALC_H