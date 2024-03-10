#if defined USE_CUDA && defined SPECTRUM

#include "cuda_integrator.h"
#include "cuda_scattering.h"
#include "global_consts.h"
#include "global_def.h"

static inline __device__ Type get_compton_frq(Type frq, Type v, Type cosf, Type cosf1, Type cosTh, Type lorenz) {

  constexpr Type hmc = kH_plank / (kM_electron * kC_Light * kC_Light);

  Type eps = hmc * frq;
  return frq * (1. - v * cosf) / (1. + eps / lorenz * (1 - cosTh) - v * cosf1);
}

static inline __device__ Type get_compton_frq(Type frq, Type cosTh) {
  constexpr Type hmc = kH_plank / (kM_electron * kC_Light * kC_Light);
  Type eps = hmc * frq;
  return frq / (1. + eps * (1.0 - cosTh));
}

static inline __device__ Type GammaC(const Type cosf, const Type v, const Vector3 &vel, const Type frq, const Type lorenz,
                                     const Vector3 &main_dir, const Vector3 &scat_dir) {

  const Type cosTh = main_dir.dot(scat_dir);
  Type frq1;
  if (v > (1. / kVelocity)) {
    const Type cosf1 = vel.dot(scat_dir) / v;
    frq1 = get_compton_frq(frq, v, cosf, cosf1, cosTh, lorenz);
  } else {
    frq1 = get_compton_frq(frq, cosTh);
  }

  constexpr Type coef = 2 * PI * kR_electron * kR_electron;
  Type frq_frac = frq1 / frq;
  Type frq1_inv = 1. / frq_frac;
  Type ds = coef * frq_frac * frq_frac * (frq_frac + frq1_inv + cosTh * cosTh - 1.0);
  return frq1_inv * ds;
}

__global__ void cuda::kernel::Get_spectrum_multi_device(const geo::grid_directions_device_t *dir, geo::grid_device_t *grid, const Type frq,
                                                        const IdType end_dir) {

  const IdType i = blockIdx.x * blockDim.x + threadIdx.x;
  const IdType k = blockIdx.y * blockDim.y + threadIdx.y;

  if (i >= grid->loc_size_gpu || k >= end_dir)
    return;

  const IdType M = dir->size;
#ifdef SINGLE_GPU
  Vector3 vel = grid->velocity[grid->shift_gpu + i]; //здесь сразу загрузили всё
#else
  Vector3 vel = grid->velocity[i];
#endif

  const Vector3 &cur_dir = dir->directions[grid->local_scattering_disp + k].dir;
  const Type *Illum = grid->illum;
  const geo::direction_device_t *all_dir = dir->directions;

  Type v2 = vel.dot(vel);

  Type lorenz = 1.0;
  Type v = 0.0;
  Type cosf = 0.0;
  if (v2 > kC_LightInv) {
    lorenz = 1. / sqrt(1. - (v2));
    v = sqrt(v2);
    cosf = vel.dot(cur_dir) / v;
  }

  Type scatter = 0;
  for (IdType num_direction = 0; num_direction < M; num_direction++) {
    Type I = Illum[i * M + num_direction];
    scatter += GammaC(cosf, v, vel, frq, lorenz, cur_dir, all_dir[num_direction].dir) * I * all_dir[num_direction].area;
  }

  ///\warning Без ослабления рассеяния нет сходимости в релятивистском режиме
  constexpr Type cfl = 0.25;

  ///\note в пределе малых частот здесь получим kSigma_thomson (размерную) (и cfl=1)
  grid->int_scattering[i * grid->local_scattering_size + k] = cfl * scatter / dir->full_area;
}

#endif //! USE_CUDA