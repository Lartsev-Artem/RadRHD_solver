/**
 * @file mpi_ext.h
 * @brief MPI расширение
 *
 * Пользовательское расширение MPI.
 * Содержит новые определённые типы данных
 */

#ifndef MPI_EXTENSION
#define MPI_EXTENSION

#include <mpi.h>

#include <cstdint>

/**
 * @brief Get the mpi id process
 *
 * @return int8_t id
 */
inline int8_t get_mpi_id();

/**
 * @brief Get the mpi count nodes
 *
 * @return int8_t count nodes
 */
inline int8_t get_mpi_np();

extern MPI_Datatype MPI_flux_t;            ///< mpi-тип для перессылки структуры ::flux_t
extern MPI_Datatype MPI_flux_illum_elem_t; ///< mpi-тип для перессылки структуры
extern MPI_Datatype MPI_hllc_value_t;      ///< mpi-тип для перессылки структуры
extern MPI_Datatype MPI_flux_all_t;        ///< mpi-тип для перессылки структуры
extern MPI_Datatype MPI_flux_elem_t;       ///< mpi-тип для перессылки структуры

#define MPI_START(argc, argv) MPI_Init(&argc, &argv);
#define MPI_END MPI_Finalize();

#endif // MPI_EXTENSION