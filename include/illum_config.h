/**
 * @file illum_config.h
 * @author Artem
 * @brief Файл конфигурации решателей излучения
 * @version 0.1
 * @date 2025-06-08
 *
 * @copyright Copyright (c) 2025
 *
 */
#if !defined ILLUM_CONFIG_H && defined ILLUM
#define ILLUM_CONFIG_H

#include "prj_config.h"

///\todo: change TRANSFER_CELL_TO_FACE and ILLUM_ON_CELL
#if 0
/* ============== Available data storage =============== */
#define ILLUM_ON_FACE 1
#define ILLUM_ON_CELL 2
/*=======================================================*/
#define ILLUM_STORAGE ILLUM_ON_FACE
#endif

#define TRANSFER_CELL_TO_FACE ///< перенос трассировки с ячеек на грани

#ifndef SPECTRUM

/// \note: энергия, поток, импульс не нужны в газодинамическом расчёте
/// для экономии памяти и времени пересылок данных можно отключить их расчёт.
/// но не для ускорения расчёта, т.к. дивергенции считаются  через эти величины
#define ON_FULL_ILLUM_ARRAYS ///< расчитывать все параметры зависящие от излучения на каждом шаге

#else
#define SAVE_FULL_SPECTRUM ///< хранить все частоты в памяти

#endif //! SPECTRUM

// #define INTERPOLATION_ON_FACES ///< расчитывать линейную 2d интерполяцию на гранях, вместо усреднения

// #define USE_TRACE_THROUGH_INNER_BOUNDARY ///< использование трассировки сквозь внутреннюю границу (доп. память и эффект луча!)

// #define ONLY_CUDA_SCATTERING ///< отключает весь расчет на видеокарте, кроме рассеяния

// #define USE_REFLECTIVE_BOUNDARY ///< включить отражающую границу

#ifdef USE_CUDA
#define SEPARATE_GPU ///< разделение данных на видеокарте

#ifdef SEPARATE_GPU

///\todo: change SINGLE_GPU and MULTI_GPU
#if 0
/* ============== Available SEPARATED modes =============== */
#define MULTI_GPU 1  ///< расчёт с применением нескольких карт на одном узле
#define SINGLE_GPU 2 ///< разделенной расчет, но с одной картой (система очереди)
/*=======================================================*/
#define SEPARATED_GPU_MODE SINGLE_GPU
#endif

#define ILLUM_ON_CELL ///< излучения хранится на ячейках, а не на гранях

// #define MULTI_GPU     ///< расчёт с применением нескольких карт на одном узле
#define SINGLE_GPU ///< разделенной расчет, но с одной картой (система очереди)

#ifdef SINGLE_GPU
#define GPU_DIV_PARAM 1 ///< параметр дробления (число <<мнимых>> карт)
#endif

#endif //! SEPARATE_GPU
#endif //! USE_CUDA

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

#if !defined SOLVERS && defined USE_CUDA
#error "bad make config. Don't use cuda without solvers!"
#endif

#endif //! ILLUM_CONFIG_H