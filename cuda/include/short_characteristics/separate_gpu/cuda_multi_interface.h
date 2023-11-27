
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

/**
 * @brief Расчёт энергии, импульса и излучения
 *
 * @param[in] id_dev номер карты (физической)
 * @param[in] im_dev номер мнимой карты
 * @param[in] grid_dir сфера направлений
 * @param[in] grid сетка
 * @param[in] st id потока
 * @return int ::e_type_completion
 */
int CalculateAllParamAsync(const int id_dev, const int im_dev, const grid_directions_t &grid_dir, grid_t &grid, e_cuda_stream_id_t st);

/**
 * @brief Расчёт интеграла рассеяния без блокировки и с асинхронной отправкой данных

 * @param[in] grid_dir сфера направлений
 * @param[in] grid сетка
 * @param[in] start_dir начало направлений для данного потока
 * @param[in] end_dir конец направлений для данного потока
 * @param[in] stream id потока
 * @return int ::e_type_completion
 */
int CalculateIntScatteringAsync(const grid_directions_t &grid_dir, grid_t &grid, const IdType start_dir, const IdType end_dir, const e_cuda_stream_id_t stream);

#ifdef MULTI_GPU
/**
 * @brief Расчёт энергии, импульса и излучения на нескольких картах
 *
 * @param[in] grid_dir сфера направлений
 * @param[in] grid сетка
 * @param[in] st id потока
 * @return int ::e_type_completion
 */
int CalculateAllParamAsync(const grid_directions_t &grid_dir, grid_t &grid, e_cuda_stream_id_t st);
#endif
} // namespace separate_device

} // namespace interface
} // namespace cuda
#endif ///!CUDA_MULTI_INTERFACE_H