/**
 * @file cuda_interface.h
 * @brief Файл хранит методы и структуры, через которые можно обращаться к через cpp код,
 * скомпилированный без nvcc
 * @note все функции модуля построены на идеалогии "В любой непонятной ситуации кладём весь процесс"
 *
 */
#if !defined CUDA_INTERFACE_H && defined USE_CUDA
#define CUDA_INTERFACE_H

#include "geo_types.h"
#include "solvers_struct.h"

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
  e_сuda_scattering_1 = 0, ///< первая часть интеграла рассеяния
  e_сuda_scattering_2 = 1, ///< вторая часть интеграла рассеяния
  e_сuda_params = 2,       ///<  расчёт физических величин (потоки, импульсы)
  // e_сuda_sender,
  e_сuda_count ///< общее число потоков
};

namespace interface {

int InitDevice(const std::string &address, const grid_directions_t &grid_dir_host, grid_t &grid_host, const int start, const int end);
void ClearDevice();
void ClearHost(grid_t &grid_host);

int CalculateAllParam(const grid_directions_t &grid_dir, grid_t &grid);
int CalculateIntScattering(const grid_directions_t &grid_dir, grid_t &grid);

void CudaWait();

} // namespace interface

} // namespace cuda

#endif //! CUDA_INTERFACE_H