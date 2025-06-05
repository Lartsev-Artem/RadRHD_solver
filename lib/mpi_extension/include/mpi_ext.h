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
extern MPI_Datatype MPI_hllc_map_t;   ///< mpi-тип для перессылки структуры ::mpi_hd_t::mpi_map_t
/**
 * @brief Структура mpi пересылок газодинамической части
 */
struct mpi_hd_t
{
  MPI_Comm comm; ///< группа на которой считается газовая часть

  struct bound_t
  {
    int64_t reg_l;    ///<  самая левая ячейка, полностью лежащая на узле
    int64_t reg_r;    ///<  самая правая ячейка, полностью лежащая на узле
    int64_t left_id;  ///< самая левая ячейка общитываемая на узле (включительно) (begin)
    int64_t right_id; ///< самая правая общитываемая на узле (не включительно) (end)
  };
  struct mpi_map_t
  {
    int64_t np_l; ///<  узел "слева"
    int64_t np_r; ///<  узел "справа"
    bound_t cell; ///< границы по ячейкам
    bound_t face; ///< границы по граням
  };
  std::vector<mpi_map_t> maps; ///< распределение ячеек по узлам

  std::vector<MPI_Request> requests_cast_phys;  ///< запросы на передачу физ. переменных
  std::vector<MPI_Request> requests_send_faces; ///< запросы на передачу потоков грани
  std::vector<MPI_Request> requests_rcv_faces;  ///< запросы на приём потоков грани
  mpi_hd_t() : comm(MPI_COMM_NULL) {}
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
