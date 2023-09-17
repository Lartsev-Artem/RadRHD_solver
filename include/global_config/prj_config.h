#ifndef PRJ_CONFIG
#define PRJ_CONFIG

/* ============== Available compiling (make) keys ===============
* CLASTER/LINUX/WINDOWS --- целевая платформа (сборка под кластер/linux/windows)
* BUILD_GRAPH          --- собирать модуль построения графов
* BUILD_DATA_TO_ILLUM  --- собирать модуль трассировки лучей
* SOLVERS              --- собирать модуль решателей
* UTILS                --- собирать модуль утилит
* USE_CUDA             --- подключение технологии cuda
* USE_MPI              --- подключение технологии mpi
*===============================================================*/

#define NUMBER_OF_MEASUREMENTS 3 // размерность решателей

/* ============== Available geometry forms =============== */
#define     Cube
#define     Step
#define     Cone
#define     Cone_JET
#define     Sphere
#define     Cylinder
/*=======================================================*/

#define GEOMETRY_TYPE Sphere


#if !defined CLASTER && !defined LINUX
#if !__NVCC__
//#define USE_VTK	  // использование vtk для вывода результатов в виде сетки
#endif
#endif

#if __NVCC__
#undef USE_MPI  // компилятор nvcc не поддерживает mpi
#endif

#if defined BUILD_GRAPH || defined BUILD_DATA_TO_ILLUM 
#define ONLY_GEO_DATA   // конфигурировать только геометричекие файлы, без трассировки
#endif


#ifdef DEBUG

#define WRITE_GLOBAL_LOG // писать лог файл

#ifdef USE_MPI
#define WRITE_MPI_LOG	// писать mpi лог файл
#endif
//#define DEBUG_MPI_RHLLC

#endif// DEBUG


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
#error "use only one platforme!!!"
#endif

#if defined WINDOWS && (defined LINUX || defined CLASTER) 
#error "use only one platforme!!!"
#endif

#if defined CLASTER && (defined WINDOWS || defined LINUX) 
#error "use only one platforme!!!"
#endif

#endif //PRJ_CONFIG
