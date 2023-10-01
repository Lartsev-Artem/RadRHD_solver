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
#include <string>

namespace cuda {

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

/**
 * @brief Безопасный вызов cuda функции c проверкой возвращаемого значения
 *
 */
#define CUDA_CALL_FUNC(_func, ...)                              \
  if (CheckError(_func(__VA_ARGS__), CONVERT_TO_STRING(_func))) \
  EXIT_ERR("Error with args: %s\n", CONVERT_TO_STRING(__VA_ARGS__))

/**
 * @brief Вызов функции ядра
 *
 */
#define CUDA_CALL_KERNEL(_func, _blocks, _threads, ...) _func<<<_blocks, _threads>>>(__VA_ARGS__)

/**
 * @brief Вызов функции ядра с указанием потока выполнения
 *
 */
#define CUDA_CALL_KERNEL_STREAM(_func, _blocks, _threads, _st, ...) _func<<<_blocks, _threads, 0, _st>>>(__VA_ARGS__)

} // namespace cuda
#endif // CUDA_DEF_H