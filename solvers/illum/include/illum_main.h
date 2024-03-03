/**
 * @file illum_main.h
 * @brief  Файл подключает модуль расчета излучения
 * @version 0.1
 * @date 2023-10-01
 *
 */
#if !defined ILLUM_MAIN_H && defined ILLUM && defined SOLVERS
#define ILLUM_MAIN_H

#include "solvers_config.h"
namespace illum {

int RunIllumModule();
int RunIllumMpiModule();

int RunIllumFacesModule();
int RunIllumMultiGpuModule();

int RunSpectrumModule(int count_states);
} // namespace illum
#endif //! ILLUM_MAIN_H