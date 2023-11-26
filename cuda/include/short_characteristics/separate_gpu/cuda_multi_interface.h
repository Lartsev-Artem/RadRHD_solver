
#if !defined CUDA_MULTI_INTERFACE_H && defined USE_CUDA
#define CUDA_MULTI_INTERFACE_H

#include "geo_types.h"
#include "solvers_struct.h"

/*! \addtogroup cuda Модуль расчёта излучения на видеокарте
    \details Модуль включает копию cpu геометрии и функции расчёта интегралов рассеяния на видеокарте.
    А также величин зависящих от излучения (потоки, импульсы)
    @{
*/

namespace cuda {
namespace interface {

/**
 * @brief Пространство имён интерфейсных функция модуля cuda для нескольких карт
 *
 */
namespace separate_device {
/**
 * @brief Инициализация видеокарты
 * \warning структура grid_host должна быть проинициализирована!
 * @param[in] grid_dir_host сфера направлений
 * @param[inout] grid_host сетка
 * @return ::e_type_completion
 */
int InitDevice(const grid_directions_t &grid_dir_host, grid_t &grid_host);

/**
 * @brief удаление структур на видеокарте
 *
 */
void ClearDevice();

/**
 * @brief удаление памяти из структуры сетки на хосте
 *
 * @param[inout] grid_host сетка
 */
void ClearHost(grid_t &grid_host);

int CalculateAllParamAsync(const grid_directions_t &grid_dir, grid_t &grid, e_cuda_stream_id_t st);

} // namespace separate_device

} // namespace interface
} // namespace cuda
#endif ///!CUDA_MULTI_INTERFACE_H