/**
 * @file cuda_ray_init.h
 * @brief Функции выделения управления памятью на видеокарте
 *
 */
#if !defined CUDA_RAY_INIT_H && defined USE_CUDA
#define CUDA_RAY_INIT_H

#include "cuda_ray_struct.h"
#include "geo_types.h"

/*! \addtogroup ray_tracing Модуль трассировки лучей, построения картинной плоскости и кривых блеска
    @{
*/

namespace cuda::ray_tracing {

/**
 * @brief Инициализация памяти на видеокарте
 *
 * @param[in] faces_host грани на хосте (с копированием)
 * @param[out] faces_device грани на видеокарте
 * @param[in] rays_host лучи на хосте(инициализированный массив. не копируется)
 * @param[out] rays_device лучи на видеокарте
 * @param[out] intersection_device коды пересечений на видеокарте
 * @return int ::e_type_completion
 */
int InitMemory(const std::vector<FaceCell> &faces_host, Face *&faces_device,
               const std::vector<Ray_t> &rays_host, Ray *&rays_device,
               int *&intersection_device);

/**
 * @brief Очистка памяти видеокарты
 *
 * @param[in] faces_device грани на видеокарте
 * @param[in] rays_device лучи на видеокарте
 * @param[in] intersection_device коды пересечений на видеокарте
 * @return int ::e_type_completion
 */
int ClearMemory(Face *&faces_device, Ray *&rays_device, int *&intersection_device);

} // namespace cuda::ray_tracing
#endif //! CUDA_RAY_INIT_H