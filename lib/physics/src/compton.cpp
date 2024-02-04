#include "global_value.h"
#include "plunk.h"

/** @brief Полное сечение рассеяния
  @param eps h*nu/(m*c^2) (для не релятивистского случая
  @todo предельные случаи
*/
static inline double sigma(double eps) {

  if (eps < 1e-4) {
    return kSigma_thomson;
  }

  if (eps > 1000) {
    constexpr double coef = PI * kR_electron * kR_electron;
    return coef / eps * (0.5 + log(2 * eps));
  }

  constexpr double coef = 2.0 * PI * kR_electron * kR_electron;
  const double eps21 = 2 * eps + 1;
  const double ln_eps21 = log(eps21) / eps;

  return coef * (((1 + eps) / (eps * eps)) * (2 * (1 + eps) / (eps21)-ln_eps21) + 0.5 * ln_eps21 - (1 + 3 * eps) / (eps21 * eps21));
}

double get_scat_coef(double frq) {
  constexpr double hmc = kH_plank / (kM_electron * kC_Light * kC_Light);
  return sigma(hmc * frq);
}

double get_scat_coef(double frq, double vel, double cosf) {
  constexpr double hmc = kH_plank / (kM_electron * kC_Light * kC_Light);
  double velC = vel / kC_Light;
  double lorenz = 1. / sqrt(1. - (velC * velC));
  return sigma(2 * hmc * frq * lorenz * (1. - velC * cosf));
}

double get_compton_frq(double frq, double vel, double cosf, double cosf1, double cosTh) {

  constexpr double hmc = kH_plank / (kM_electron * kC_Light * kC_Light);
  double velC = vel / kC_Light;
  double lorenz = 1. / sqrt(1. - (velC * velC));
  double eps = hmc * frq;

  return frq * (1. - velC * cosf) / (1. + eps / lorenz * (1 - cosTh) - velC * cosf1);
}

double get_dif_scat_coef(double frq, double vel, double cosf, double cosf1, double cosTh) {
  constexpr double coef = 0.5 * kR_electron * kR_electron;
  double frq1 = get_compton_frq(frq, vel, cosf, cosf1, cosTh) / frq;
  double frq1_inv = 1. / frq1;
  return coef * frq1 * frq1 * (frq1 + frq1_inv + cosTh * cosTh - 1);
}

Type get_int_func(const Type frq, const Vector3 &vel, const Vector3 &dir, const Vector3 &dir_scat, const Type *Illum) {

  const double v = vel.norm();
  const double cosTh = dir.dot(dir_scat);
  const double cosf1 = vel.dot(dir_scat) / v;
  const double cosf = vel.dot(dir) / v;

  Type frq1 = get_compton_frq(frq, v, cosf, cosf1, cosTh);

  // Type ds = get_dif_scat_coef(frq, v, cosf, cosf1, cosTh);
  constexpr double coef = 0.5 * kR_electron * kR_electron;
  double frq_frac = frq1 / frq;
  double frq1_inv = 1. / frq_frac;
  Type ds = coef * frq_frac * frq_frac * (frq_frac + frq1_inv + cosTh * cosTh - 1.0);

  int idx = get_frq_idx(frq1);

  return frq1_inv * ds * Illum[idx];
}

#if 0 // test for Wolfram plot
#include <fstream>
#include <iomanip>
using namespace std;
int main() {

  constexpr double kSigma_thomson = 6.65210 * 1e-29;

  ofstream ofile;
  ofile.open("compton.txt");

  const double frq = kM_electron * kC_Light * kC_Light / kH_plank;

  const Vector3 vel(0.99 * kC_Light, 0, 0);
  const Vector3 dir(1, 0, 0);
  const Vector3 dir_scat(0, 1, 0);

  // const double v = vel.norm();
  // const double cosTh = dir.dot(dir_scat);
  // const double cosf1 = vel.dot(dir_scat) / v;
  // const double cosf = vel.dot(dir) / v;

  const double v = 0.7 * kC_Light;
  double cosTh = 0.7;
  double cosf1 = 0.0;
  double cosf = cos(0.);

  // /// полное сечение
	// std::vector<Type> freq;
	// get_splitting_spectrum(freq);   

  //  for (int i=freq.size()-1; i>=0 ; i--) {
	// cosf=0;
  //   ofile << setprecision(16) << freq[i] << ' '
  //         << get_scat_coef(freq[i],v,cosf)/kSigma_thomson
  //         << '\n';
  // }

  for (double i = 0; i < PI; i += (PI / 360)) {
    cosTh = cos(i);
    cosf1 = cos(i - 0); // Cos[\[Theta] - \[Phi]]
    ofile << setprecision(16) << i << ' '
          << (get_dif_scat_coef(frq, v, cosf, cosf1, cosTh) / (get_dif_scat_coef(frq, v, cosf, cos(0.), cos(0.))))
          << '\n';
  }
  ofile.close();

  return 0;
}
#endif