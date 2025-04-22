/**
 * @file mpi_ext.h
 * @brief MPI расширение
 *
 * Пользовательское расширение MPI.
 * Содержит новые определённые типы данных
 */

#ifndef MPI_EXTENSION
#define MPI_EXTENSION

#include <cstdint>

/**
 * @brief Get the mpi id process
 *
 * @return int8_t id
 */
/*inline*/ int8_t get_mpi_id();

/**
 * @brief Get the mpi count nodes
 *
 * @return int8_t count nodes
 */
/*inline*/ int8_t get_mpi_np();

#ifdef USE_MPI
#include <global_types.h>
#include <mpi.h>
#include <vector>
/**
 * @brief Get the mpi id process on comm
 *
 * @param[in] comm - communicator
 * @return int8_t id
 */
int8_t get_mpi_id(const MPI_Comm &comm);

/**
 * @brief Get the mpi count nodes on comm
 *
 * @param[in] comm - communicator
 * @return int8_t count nodes
 */
int8_t get_mpi_np(const MPI_Comm &comm);

extern MPI_Datatype MPI_phys_val_t;   ///< mpi-тип для перессылки структуры ::elem_t::phys_val
extern MPI_Datatype MPI_flux_t;       ///< mpi-тип для перессылки структуры ::flux_t
extern MPI_Datatype MPI_hllc_value_t; ///< mpi-тип для перессылки структуры ::hllc_value_t (динамический расчёт)
extern MPI_Datatype MPI_flux_elem_t;  ///< mpi-тип для перессылки структуры ::elem_t::phys_val ,elem_t::conv_val

/**
 * @brief Структура mpi пересылок газодинамической части
 */
struct mpi_hllc_t {
  MPI_Comm comm;                         ///< группа на которой считается газовая часть
  std::vector<IdType> send_cells;        ///< кол-во отправок для каждого процесса
  std::vector<IdType> disp_cells;        ///< смещения по ячейкам для процессов
  std::vector<IntId> id_irregular_faces; ///< номера граней с границе областей процессов

  std::vector<MPI_Request> requests_cast_phys;  ///< запросы на передачу физ. переменных
  std::vector<MPI_Request> requests_send_faces; ///< запросы на передачу потоков грани
  std::vector<MPI_Request> requests_rcv_faces;  ///< запросы на приём потоков грани
  mpi_hllc_t() : comm(MPI_COMM_NULL) {}
};

/**
 * @brief Инициализация структур в MPI_TYPE
 *
 */
void InitMPiStruct();

#define MPI_START(argc, argv) MPI_Init(&argc, &argv);
#define MPI_END MPI_Finalize();

#define MPI_BARRIER(comm) MPI_Barrier(comm);

#else
#define CONVERT_TO_STRING(s, ...) #s #__VA_ARGS__
#define MPI_START(argc, argv) CONVERT_TO_STRING(argc, argv);
#define MPI_END
#define MPI_BARRIER(comm) CONVERT_TO_STRING(comm);
#endif
#endif // MPI_EXTENSION
