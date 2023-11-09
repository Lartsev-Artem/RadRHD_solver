#if !defined ILLUM_MPI_STRUCT_H && defined ILLUM && defined SOLVERS && defined USE_MPI
#define ILLUM_MPI_STRUCT_H

#include "geo_types.h"

extern MPI_Datatype MPI_RECV_ILLUM_T; ///< структура дл отправки по направлению, а приему по ячейкам

namespace illum {

void MpiInitStruct(const grid_directions_t &grid);
} // namespace illum
#endif //! ILLUM_MPI_STRUCT_H