/**
 * @file solvers_config.h
 * @brief Файл конфигурации проекта
 *
 * \details Доступны следующие конфигурации:
 * \arg HLLC газодинамический решатель на основе метода Годунова с аппроксимацией потоков методом hllc
 * \arg RHLLC релятивисткий газодинамический решатель. Обобщение hllc
 * \arg ILLUM модуль расчёта излучения
 * \note в ручном режиме можно отключить макросы ILLUM_MPI, HLLC_MPI, RHLLC_MPI
 * для генерации отдельных решателей без поддержки MPI.
 */

#if !defined SOLVE_CONFIG_H && defined SOLVERS
#define SOLVE_CONFIG_H

#include "prj_config.h"

//#define HLLC
//#define RHLLC
// #define ILLUM

/// \note: энергия, поток, импульс не нужны в газодинамическом расчёте
/// для экономии памяти и времени пересылок данных можно отключить их расчёт.
/// но не для ускорения расчёта, т.к. дивергенции считаются  через эти величины
#ifdef ILLUM
#define ON_FULL_ILLUM_ARRAYS ///< расчитывать все параметры зависящие от излучения на каждом шаге
// #define INTERPOLATION_ON_FACES ///< расчитывать линейную 2d интерполяцию на гранях, вместо усреднения
#endif

#ifdef USE_MPI

#ifdef ILLUM
#define ILLUM_MPI ///< расчёт излучения с технологией MPI
#endif

#ifdef HLLC
#define HLLC_MPI ///< расчёт hllc с технологией MPI
#endif

#ifdef RHLLC
#define RHLLC_MPI ///< расчёт реляивиской постановки с технологией MPI
#endif

#endif //! USE_MPI

/* ============== Check config project =============== */

#if defined USE_CUDA && !defined ILLUM
#error "Bad config. CUDA only with ILLUM!!!"
#endif

#if NUMBER_OF_MEASUREMENTS != 3 && defined ILLUM
#error "Bad config. ILLUM is only available in 3D"
#endif

#if defined HLLC && defined RHLLC
#error "Bad config. There can be only one task at a time (HLLC or RHLLC)"
#endif

#endif //! SOLVE_CONFIG_H