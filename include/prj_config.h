/**
 * @file prj_config.h
 * @brief Файл конфигурации проекта
 * @details
 * Available compiling (make) keys:
 * @arg BUILD_GRAPH           --- собирать модуль построения графов
 * @arg MAKE_TRACE            --- собирать модуль трассировки лучей
 * @arg SOLVERS               --- собирать модуль решателей
 * @arg UTILS                 --- собирать модуль утилит
 * @arg USE_CUDA              --- подключение технологии cuda
 * @arg USE_MPI               --- подключение технологии mpi
 * @arg RAD_RHD               --- радиационная релятивистская газовая динамика
 * @arg SPECTRUM              --- построение спектра путем решения уравнения переноса излучения
 * @arg DEBUG                 --- сборка с логами
 *===============================================================* /
 * @version 0.1
 * @date 2023-09-21
 *
 * @copyright Copyright (c) 2023
 *
 */

/** namespace Rrhd - Relativistic radiation hydrodynamic  */

#ifndef PRJ_CONFIG
#define PRJ_CONFIG

#define NUMBER_OF_MEASUREMENTS 3 ///< размерность решателей

#if !__NVCC__
// #define USE_VTK ///< использование vtk для вывода результатов
#else
#undef USE_MPI //! компилятор nvcc не поддерживает mpi
#endif

#ifdef DEBUG

// #define LOG_OUT_TO_SCREEN ///<выводить лог не в файл, а на экран

#define WRITE_GLOBAL_LOG ///< писать лог файл

#ifdef USE_MPI
#define WRITE_MPI_LOG // писать mpi лог файл
#endif

#endif // DEBUG

#endif // PRJ_CONFIG
