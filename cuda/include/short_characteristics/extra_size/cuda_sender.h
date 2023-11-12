#if !defined CUDA_ILLUM_SENDER_H && defined USE_CUDA
#define CUDA_ILLUM_SENDER_H

#include "cuda_struct.h"
#include "solvers_struct.h"

void InitGridOnDeviceMultiCuda(int id_dev, const grid_t &grid_host,
                               const IdType size_dir, const IdType start_dir, const IdType end_dir,
                               cuda::geo::grid_device_t *&grid_device);

#endif //! CUDA_ILLUM_SENDER_H