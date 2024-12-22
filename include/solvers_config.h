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

#define TRANSFER_CELL_TO_FACE ///< перенос трассировки с ячеек на грани

/// \note: энергия, поток, импульс не нужны в газодинамическом расчёте
/// для экономии памяти и времени пересылок данных можно отключить их расчёт.
/// но не для ускорения расчёта, т.к. дивергенции считаются  через эти величины
#ifdef ILLUM

#ifndef SPECTRUM
#define ON_FULL_ILLUM_ARRAYS ///< расчитывать все параметры зависящие от излучения на каждом шаге
#endif
//#define INTERPOLATION_ON_FACES ///< расчитывать линейную 2d интерполяцию на гранях, вместо усреднения

// #define USE_TRACE_THROUGH_INNER_BOUNDARY ///< использование трассировки сквозь внутреннюю границу (доп. память и эффект луча!)

// #define ONLY_CUDA_SCATTERING ///< отключает весь расчет на видеокарте, кроме рассеяния

//#define USE_REFLECTIVE_BOUNDARY ///< включить отражающую границу

#ifdef SPECTRUM
#define SAVE_FULL_SPECTRUM ///< хранить все частоты в памяти
#endif

#ifdef USE_CUDA
#define SEPARATE_GPU ///< разделение данных на видеокарте
#endif

//#define SPECTRUM ///< перенос в полной постановке, без серой материи

#ifdef SEPARATE_GPU
#define ILLUM_ON_CELL ///< излучения хранится на ячейках, а не на гранях
// #define MULTI_GPU     ///< расчёт с применением нескольких карт на одном узле
#define SINGLE_GPU ///< разделенной расчет, но с одной картой (система очереди)

#ifdef SINGLE_GPU
#define GPU_DIV_PARAM 1 ///< параметр дробления (число <<мнимых>> карт)
#endif
#endif

#endif //! ILLUM

#ifdef USE_MPI

#ifdef ILLUM
#define ILLUM_MPI ///< расчёт излучения с технологией MPI
#endif

#if (defined HLLC) || (defined RHLLC)

#define IDEAL 1 ///< Eos: Г-law
#define TAUB 2  ///< Eos: Taub-Mathews

#define EOS TAUB

#endif

#endif //! USE_MPI

/* ============== Check config project =============== */

#if defined MULTI_GPU && defined SINGLE_GPU
#error "Bad config. A single type of separation is needed !!!"
#endif

#if defined USE_CUDA && !defined ILLUM
#error "Bad config. CUDA only with ILLUM!!!"
#endif

#if NUMBER_OF_MEASUREMENTS != 3 && defined ILLUM
#error "Bad config. ILLUM is only available in 3D"
#endif

#if defined HLLC && defined RHLLC
#error "Bad config. There can be only one task at a time (HLLC or RHLLC)"
#endif

#if defined RAD_RHD && defined SPECTRUM
#error "Bad config. There can be only one task at a time (SPECTRUM or RAD_RHD)"
#endif

#if defined RAD_RHD && (!defined RHLLC || !defined ILLUM)
#error "Bad config. RAD_RHD is only available with RHLLC and ILLUM"
#endif

#endif //! SOLVE_CONFIG_H