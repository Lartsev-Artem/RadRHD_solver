#ifndef MPI_EXTENSION
#define MPI_EXTENSION

#include <stdint.h>
#include <mpi.h>

inline int8_t get_mpi_id();
inline int8_t get_mpi_np();

extern MPI_Datatype MPI_flux_t;
extern MPI_Datatype MPI_flux_illum_elem_t;
extern MPI_Datatype MPI_hllc_value_t;
extern MPI_Datatype MPI_flux_all_t;
extern MPI_Datatype MPI_flux_elem_t;

#define MPI_START(argc, argv) MPI_Init(&argc, &argv);
#define MPI_END MPI_Finalize();

#endif //MPI_EXTENSION