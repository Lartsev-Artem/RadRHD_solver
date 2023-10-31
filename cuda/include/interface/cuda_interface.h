/**
 * @file cuda_interface.h
 * @brief Файл хранит методы и структуры, через которые можно обращаться через cpp код,
 * скомпилированный без nvcc
 * @note все функции модуля построены на идеалогии "В любой непонятной ситуации кладём весь процесс"
 *
 */
#if !defined CUDA_INTERFACE_H && defined USE_CUDA
#define CUDA_INTERFACE_H

#include "geo_types.h"
#include "solvers_struct.h"

/*! \addtogroup cuda Модуль расчёта излучения на видеокарте
    \details Модуль включает копию cpu геометрии и функции расчёта интегралов рассеяния на видеокарте.
    А также величин зависящих от излучения (потоки, импульсы)
    @{
*/

/**
 * @brief Пространство имён модуля расчета на видеокарте
 *
 */
namespace cuda {

/**
 * @brief код cuda потока
 *
 */
enum e_cuda_stream_id_t {
  e_cuda_scattering_1 = 0, ///< первая часть интеграла рассеяния
  e_cuda_scattering_2 = 1, ///< вторая часть интеграла рассеяния
  e_cuda_params = 2,       ///<  расчёт физических величин (потоки, импульсы)
  // e_cuda_sender,
  e_cuda_streams_count ///< общее число потоков
};

/**
 * @brief Пространство имён интерфейсных функция модуля cuda
 *
 */
namespace interface {

/**
 * @brief Инициализация видеокарты
 * \details установка номера видеокарты, инициализация сеток и выделения  host памяти в cpu сетки
 * \warning структура grid_host должна быть проинициализирована!
 * @param[in] address путь до файлов геометрии
 * @param[in] grid_dir_host сфера направлений
 * @param[inout] grid_host сетка
 * @return ::e_type_completion
 */
int InitDevice(const std::string &address, const grid_directions_t &grid_dir_host, grid_t &grid_host);

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
 * @brief Расчёт энергии, импульса, потока излучения
 * @note с блокирующем копированием дивергенций и асинхронном копированием
 * потоков, энергий, импульсов
 *
 * @param[in] grid_dir  сфера направлений
 * @param[out] grid сетка
 * @return ::e_type_completion
 */
int CalculateAllParam(const grid_directions_t &grid_dir, grid_t &grid);

/**
 * @brief Расчёт энергии, импульса, потока излучения
 *
 * @param[in] grid_dir сфера направлений
 * @param[inout] grid сетка
 * @param[in] st id потока
 * @return int ::e_type_completion
 * @warning Перед чтением данных на хосте необходимо синхронизироваться с потоком
 */
int CalculateAllParamAsync(const grid_directions_t &grid_dir, grid_t &grid, e_cuda_stream_id_t st);
/**
 * @brief Расчёт интеграла рассеяния (с копированием излучения на карту)
 *
 * @param[in] grid_dir  сфера направлений
 * @param[out] grid сетка
 * @return ::e_type_completion
 */
int CalculateIntScattering(const grid_directions_t &grid_dir, grid_t &grid);

/**
 * @brief Расчёт интеграла рассеяния без блокировки и с асинхронной отправкой данных
 *
 * @param[in] grid_dir сфера направлений
 * @param[in] grid сетка
 * @param[in] start_dir начало направлений для данного потока
 * @param[in] end_dir конец направлений для данного потока
 * @param[in] stream id потока
 * @return int ::e_type_completion
 */
int CalculateIntScatteringAsync(const grid_directions_t &grid_dir, grid_t &grid, const IdType start_dir, const IdType end_dir, const e_cuda_stream_id_t stream);

/**
 * @brief Барьерная синхронизация с cuda
 *
 */
void CudaWait();

/**
 * @brief Синхронизация с отдельным потоком
 *
 * @param[in] stream_id id потока
 */
void CudaSyncStream(const e_cuda_stream_id_t stream_id);

/**
 * @brief Инициализация потоков cuda с приоритетами
 *
 */
void SetStreams();

/**
 * @brief Не блокирующая отправка данных на видеокарту
 *
 * @param[in] size размер данных (в элементах)
 * @param[in] shift сдвиг в массиве
 * @param[in] Illum_host массив излучения
 */
void CudaSendIllumAsync(const IdType size, const IdType shift, const Type *Illum_host);

} // namespace interface

} // namespace cuda

#endif //! CUDA_INTERFACE_H