/**
 * @file illum_mpi_sender.h
 * @brief Вспомогательные структуры для mpi обмена
 *
 */
#if !defined ILLUM_MPI_SENDER_H && defined ILLUM && defined SOLVERS && defined USE_MPI
#define ILLUM_MPI_SENDER_H

#include "geo_types.h"
#include "solvers_struct.h"

/*! \addtogroup illum Модуль расчёта излучения
    @{
*/

namespace illum {

/// @brief структура для организации неблокирующих пересылок
struct mpi_sender_t {
  IdType size;                            ///< размер секции по направлениям (не по запросам)
  std::vector<MPI_Request> requests_rcv;  //все запросы сообщений отправки и принятия
  std::vector<MPI_Status> status_rcv;     //статусы всех обменов
  std::vector<int> flags_send_to_gpu;     //флаги указывающие на отправку пакета на gpu
  std::vector<MPI_Request> requests_send; //все запросы сообщений отправки и принятия
};

namespace separate_gpu {
/**
 * @brief Инициализация MPI обмена, для типа хранения по ячейкам с раздельным расчетом на видеокарте
 *
 * @param[in] comm коммуниктор для излучения (для случая независимого счета rhd и излучения)
 * @param[in] grid_dir сетка направлений
 * @param[in] grid сетка пространственная
 */
void InitSender(const MPI_Comm &comm, const grid_directions_t &grid_dir, const grid_t &grid);
} // namespace separate_gpu

namespace gpu_async {

/**
 * @brief Инициализация MPI обмена, для типа хранения по направлениям с полным хранением излучения на карте
 *
 * @param[in] grid_dir сетка направлений
 * @param[in] grid сетка пространственная
 */
void InitSender(const grid_directions_t &grid_dir, const grid_t &grid);
} // namespace gpu_async
} // namespace illum

extern illum::mpi_sender_t section_1;
extern illum::mpi_sender_t section_2;
extern MPI_Comm MPI_COMM_ILLUM;
extern std::vector<IdType> disp_illum;

#endif //! ILLUM_MPI_SENDER_H