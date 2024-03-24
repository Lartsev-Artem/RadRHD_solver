/**
 * @file cuda_def.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2023-10-01
 *
 * @copyright Copyright (c) 2023
 *
 */
#if !defined CUDA_DEF_H && defined USE_CUDA
#define CUDA_DEF_H

#include "dbgdef.h"
#include "solvers_config.h"
#include <string>

/*! \addtogroup cuda Модуль расчёта излучения на видеокарте
    @{
*/

namespace cuda {

#define CONVERT_TO_STRING(s, ...) #s #__VA_ARGS__

/**
 * @brief Функция проверки работы функций cuda с расшифровкой кода ошибки
 *
 * @tparam Err_t
 * @param cudaStatus код ошибки
 * @param my_text доп. текст сообщения
 * @return ::e_type_completion
 */
template <typename Err_t>
int inline CheckError(Err_t cudaStatus, const std::string my_text = "") {
  if (cudaStatus != cudaSuccess) {
    WRITE_LOG_ERR("Error cuda code: %d \n%s \n%s\n", cudaStatus, cudaGetErrorString(cudaStatus), my_text.c_str());
    return e_completion_fail;
  }
  return e_completion_success;
}

#ifdef DEBUG
/**
 * @brief Безопасный вызов cuda функции c проверкой возвращаемого значения
 *
 */
#define CUDA_CALL_FUNC(_func, ...)                              \
  if (CheckError(_func(__VA_ARGS__), CONVERT_TO_STRING(_func))) \
  EXIT_ERR("Error with args: %s\n", CONVERT_TO_STRING(__VA_ARGS__))
#else
#define CUDA_CALL_FUNC(_func, ...) _func(__VA_ARGS__);
#endif

#ifdef SPECTRUM
#ifdef SAVE_FULL_SPECTRUM
#define BS 4 ///< размер блока
#else
#define BS 16 /// \note на nvidia gtx 1080 6gb почему не хватает ресурсов для запуска
#endif
#else
#define BS 32 ///< размер блока
#endif

#define CUDA_BLOCKS_3D(val, cells, dirs, frqs) dim3 val((cells + BS - 1) / BS, (dirs + BS - 1) / BS, (frqs + BS - 1) / BS);
#define CUDA_TREADS_3D(val) dim3 val(BS, BS, BS);

#define CUDA_BLOCKS_2D(val, cells, dirs) dim3 val((cells + BS - 1) / BS, (dirs + BS - 1) / BS);
#define CUDA_TREADS_2D(val) dim3 val(BS, BS);

#define CUDA_BLOCKS_1D(val, cells) dim3 val((cells + BS - 1) / BS);
#define CUDA_TREADS_1D(val) dim3 val(BS);

} // namespace cuda
#endif // CUDA_DEF_H