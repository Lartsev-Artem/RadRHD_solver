/**
 * @file rhllc_3d_mpi.h
 * @author Artem
 * @brief rhllc solver with mpi support
 * @version 0.1
 * @date 2025-06-10
 *
 * @copyright Copyright (c) 2025
 *
 */
#if !defined RHLLC_MPI_H && defined RHLLC && defined USE_MPI
#define RHLLC_MPI_H

#include "mpi_ext.h"
#include "solvers_struct.h"

/*! \addtogroup rhllc Модуль расчета газовой динамики в релятивистской постановке
    @{
*/

namespace rrhd {
namespace rhllc_mpi {

void Hllc3d(const Type tau, grid_t &grid);

}; // namespace rhllc_mpi

}; // namespace rrhd

#endif //! RHLLC_MPI_H