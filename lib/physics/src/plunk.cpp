#include "plunk.h"

#include "global_value.h"

namespace ksplit {
constexpr double betta = 35000;
constexpr double frq0 = 1e24;
constexpr int N = 1000;
} // namespace ksplit

void get_splitting_spectrum(std::vector<Type> &spectrum) {

  spectrum.resize(ksplit::N);
  for (int x = 0; x < ksplit::N; x++) {
    spectrum[x] = ksplit::frq0 * exp(-x * x / ksplit::betta);
  }
}

int get_frq_idx(double frq) {
  if (frq > ksplit::frq0) {
    return ksplit::N - 1;
  }
  return floor(sqrt(-ksplit::betta * log(frq / ksplit::frq0)));
}

static inline double q1(double x) {
  return (x * x * x + 3 * x * x + 6 * x + 7.28);
}
static inline double q0(double x) {
  return (x * x * x * (1. / 3. - x / 8 + x * x / 62.4));
}

/**
 * @brief Приближенное вычисление интеграла Планка предложенное Гольдиным (ost Четверушкин)
 */
static inline double f_Goldin(double x) {
  if (x <= 2.0) {
    return q0(x);
  }

  return (PI * PI * PI * PI) / 15 - exp(-x) * q1(x);
}

static double f_Goldin_sigma_log(double T, double nu, double nu0) {

  if (nu0 > nu) {
    abort();
  }

  constexpr double k_plankB = 2 * (k_boltzmann * k_boltzmann * k_boltzmann * k_boltzmann) / (kC_Light * kC_Light * kH_plank * kH_plank * kH_plank);
  constexpr double k_Log_plankB = -12.7917563178452838;
  constexpr double c = kH_plank / k_boltzmann;
  const double coef = c / T;
  const double x1 = coef * nu;
  const double x0 = coef * nu0;

  if (x1 <= 2.0) // && x0
  {
    return k_Log_plankB + log(q0(x1) - q0(x0));
  } else if (x1 > 2.0 && x0 < 2.0) {
    return k_Log_plankB + log(((PI * PI * PI * PI) / 15 - q0(x0)) - exp(-x1) * q1(x1)); //основной вклад от x0
  } else                                                                                // x1>2 && x0>2
  {

    const double y = 0.5 * coef * (nu - nu0);

    /// \warning реально переход в первое условие намного раньше т.к. операция [exp(x)-exp(-x)] без потери разрядов
    /// возможно только при |x|<37 (но тут exp(x)*q(x))

    if (y > 500) // граница перехода exp(x) в inf. (реально y>700. С запасом)
    {
      return k_Log_plankB - 0.5 * coef * (nu + nu0) + y + log(q1(x0)); //отбросили вычитаемое, взяли логарифм полуаналитически
    } else {
      return k_Log_plankB - 0.5 * coef * (nu + nu0) + log(exp(y) * q1(x0) - exp(-y) * q1(x1));
    }
  }
}

double B_Plank(double T, double nu, double nu0) {
  return (T * T * T * T) * exp(f_Goldin_sigma_log(T, nu, nu0));
}

double B_Plank_log(double T, double nu, double nu0) {
  return 4 * log(T) + (f_Goldin_sigma_log(T, nu, nu0));
}

#if 0 // test for Wolfram plot
#include <fstream>
#include <iomanip>
using namespace std;
int main() {

  ofstream ofile;
  ofile.open("plank.txt");

  double T = 270;
  double d_nu = 1e14;
  double nu0 = 3 * 1e14; // 1000нм
  double nu = 1e17;      //

  std::vector<Type> frq;
  get_splitting_spectrum(frq);

  for (size_t i = 0; i < frq.size() - 1; i++) {
    ofile << setprecision(16) << frq[i] / frq.back() << ' ' << (B_Plank_log(T, frq[i], frq[i + 1])) << '\n';
  }

  ofile << setprecision(16) << 0 << ' ' << (B_Plank_log(T, frq.back(), 0)) << '\n';
  ofile.close();
  return 0;
}
#endif