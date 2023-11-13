#if !defined ILLUM_MPI_SENDER_H && defined ILLUM && defined SOLVERS && defined USE_MPI
#define ILLUM_MPI_SENDER_H

#include "geo_types.h"
#include "solvers_struct.h"

namespace illum {

struct mpi_sender_t {
  IdType size;                            ///< размер секции по направлениям (не по запросам)
  std::vector<MPI_Request> requests_rcv;  //все запросы сообщений отправки и принятия
  std::vector<MPI_Status> status_rcv;     //статусы всех обменов
  std::vector<int> flags_send_to_gpu;     //флаги указывающие на отправку пакета на gpu
  std::vector<MPI_Request> requests_send; //все запросы сообщений отправки и принятия
};

namespace extra_size {
void InitSender(const MPI_Comm &comm, const grid_directions_t &grid_dir, const grid_t &grid);
}

namespace gpu_async {
void InitSender(const grid_directions_t &grid_dir, const grid_t &grid);
}
} // namespace illum

extern illum::mpi_sender_t section_1;
extern illum::mpi_sender_t section_2;
extern MPI_Comm MPI_COMM_ILLUM;
extern std::vector<IdType> disp_illum;

#endif //! ILLUM_MPI_SENDER_H