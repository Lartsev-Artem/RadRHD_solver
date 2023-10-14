#ifdef USE_CUDA
#include "cuda_def.h"
#include "cuda_interface.h"

#include "cuda_struct.h"

void cuda::interface::CudaWait() {
  CUDA_CALL_FUNC(cudaDeviceSynchronize);
}

cudaStream_t cuda_streams[cuda::e_сuda_streams_count];

void cuda::interface::CudaSyncStream(const e_cuda_stream_id_t stream_id) {
  CUDA_CALL_FUNC(cudaStreamSynchronize, cuda_streams[stream_id]);
}

void cuda::interface::SetStreams() {

  // get the range of stream priorities for this device
  int priority_high, priority_low;
  cudaDeviceGetStreamPriorityRange(&priority_low, &priority_high);

  CUDA_CALL_FUNC(cudaStreamCreateWithPriority, &cuda_streams[e_сuda_scattering_1], cudaStreamNonBlocking, priority_high);
  CUDA_CALL_FUNC(cudaStreamCreateWithPriority, &cuda_streams[e_сuda_scattering_2], cudaStreamNonBlocking, (priority_high + priority_low) / 2);
  CUDA_CALL_FUNC(cudaStreamCreateWithPriority, &cuda_streams[e_сuda_params], cudaStreamNonBlocking, priority_low);
}

void cuda::interface::CudaSendIllumAsync(const int size, const int shift, const Type *Illum_host) {

  CUDA_CALL_FUNC(cudaMemcpyAsync, device_host_ptr.illum + shift, Illum_host + shift, size * sizeof(Illum_host[0]),
                 cudaMemcpyHostToDevice, cuda_streams[e_сuda_scattering_1]);
}

#endif //! USE_CUDA

#if 0
static cudaEvent_t cuda_events[cuda::e_сuda_streams_count * 3];
void CudaEventPrint() {
  for (int i = 0; i < cuda::e_сuda_streams_count * 3; i++) {
    // synchronize
    cudaEventSynchronize(cuda_events[i]); // optional
  }

  float dt_ms;
  cudaEventElapsedTime(&dt_ms, cuda_events[0], cuda_events[1]);
  printf("time kernel1 %f\n", dt_ms);

  cudaEventElapsedTime(&dt_ms, cuda_events[1], cuda_events[2]);
  printf("time send1 %f\n", dt_ms);

  cudaEventElapsedTime(&dt_ms, cuda_events[3], cuda_events[4]);
  printf("time kernel2 %f\n", dt_ms);

  cudaEventElapsedTime(&dt_ms, cuda_events[4], cuda_events[5]);
  printf("time send2 %f\n\n", dt_ms);
}
#endif