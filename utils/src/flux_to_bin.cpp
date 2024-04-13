#include "utils.h"

#include "global_types.h"
#include "global_value.h"
#include "solvers_struct.h"

#include "reader_bin.h"
#include "writer_bin.h"

int FUNC_NAME(FluxToBin)(int argc, char *argv[]) {
  if (argc < 3) {
    printf("Error input data!\n");
    printf("Input:\n");
    printf("address_solve (like \"path\\file\")\n");
    printf("max_number_of_iter\n");
    printf("Add params: remove exist files [0/1]\n");
    return e_completion_fail;
  }

  const std::string address_solve = argv[1];
  const int max_number_of_iter = std::stoi(argv[2]);

  bool del = false;
  if (argc == 4) {
    del = std::stoi(argv[3]);
  }

  std::vector<flux_t> flux;
  std::vector<Type> den;
  std::vector<Type> pres;
  std::vector<Vector3> vel;
  for (int i = 0; i < max_number_of_iter; i++) {
    if (files_sys::bin::ReadSimple(address_solve + std::to_string(i) + F_FLUX, flux)) {
      return e_completion_fail;
    }

    den.resize(flux.size());
    pres.resize(flux.size());
    vel.resize(flux.size());

    for (size_t j = 0; j < flux.size(); j++) {
      den[j] = flux[j].d;
      pres[j] = flux[j].p;
      vel[j] = flux[j].v;
    }

    files_sys::bin::WriteSimple(address_solve + std::to_string(i) + F_DENSITY, den);
    files_sys::bin::WriteSimple(address_solve + std::to_string(i) + F_PRESSURE, pres);
    files_sys::bin::WriteSimple(address_solve + std::to_string(i) + F_VELOCITY, vel);

    if (del) {
      std::remove((address_solve + std::to_string(i) + F_FLUX).c_str());
    }
  }
  return e_completion_success;
}