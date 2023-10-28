/**
 * @file prj_config.h
 * @brief Файл конфигурации проекта
 * @details
 * Available compiling (make) keys:
 * @arg CLASTER/LINUX/WINDOWS --- целевая платформа (сборка под кластер/linux/windows)
 * @arg BUILD_GRAPH           --- собирать модуль построения графов
 * @arg MAKE_TRACE            --- собирать модуль трассировки лучей
 * @arg SOLVERS               --- собирать модуль решателей
 * @arg UTILS                 --- собирать модуль утилит
 * @arg USE_CUDA              --- подключение технологии cuda
 * @arg USE_MPI               --- подключение технологии mpi
 * @arg DEBUG                 --- сборка с логами
 *===============================================================* /
 * @version 0.1
 * @date 2023-09-21
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef PRJ_CONFIG
#define PRJ_CONFIG

#define NUMBER_OF_MEASUREMENTS 3 ///< размерность решателей

/* ============== Available geometry forms =============== */
#define Cube 1
#define Step 2
#define Cone 3
#define Cone_JET 4
#define Sphere 5
#define Cylinder 6
#define TEST_ELLIPSE 7
#define MAIN_ELLIPSE 8
/*=======================================================*/

#define GEOMETRY_TYPE Cube // MAIN_ELLIPSE // Sphere

#if !defined CLASTER
#if !__NVCC__
//#define USE_VTK ///< использование vtk для вывода результатов в виде
//сетки
#endif
#endif

#if __NVCC__
#undef USE_MPI //! компилятор nvcc не поддерживает mpi
#endif

#if defined BUILD_GRAPH || defined BUILD_DATA_TO_ILLUM
//#define ONLY_GEO_DATA ///< конфигурировать только геометричеcкие файлы, без трассировки
#endif

#ifdef DEBUG

//#define LOG_OUT_TO_SCREEN ///<выводить лог не в файл, а на экран

#define WRITE_GLOBAL_LOG ///< писать лог файл

#ifdef USE_MPI
#define WRITE_MPI_LOG // писать mpi лог файл
#endif
//#define DEBUG_MPI_RHLLC

#endif // DEBUG

/// \todo: delete it!!!
// #if defined UTILS && !defined SOLVERS
// #define ILLUM
// #define HLLC
// #endif //UTILS

/* ============== Check config project =============== */

#if !defined SOLVERS && defined USE_CUDA
#error "bad make config. Don't use cuda without solvers!"
#endif

#if defined LINUX && (defined WINDOWS || defined CLASTER)
#error "use only one platform!!!"
#endif

#if defined WINDOWS && (defined LINUX || defined CLASTER)
#error "use only one platform!!!"
#endif

#if defined CLASTER && (defined WINDOWS || defined LINUX)
#error "use only one platform!!!"
#endif

#endif // PRJ_CONFIG
