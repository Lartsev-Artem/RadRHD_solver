/**
 * @file rhllc_mpi.h
 * @version 0.1
 * @date 2024-02-15
 *
 */
#if !defined RHLLC_MPI_H && defined SOLVERS && defined USE_MPI
#define RHLLC_MPI_H

#include "solvers_struct.h"

namespace rhllc_mpi {
void Hllc3dStab(const Type tau, grid_t &grid);
} // namespace rhllc_mpi

#endif //! RHLLC_MPI_H