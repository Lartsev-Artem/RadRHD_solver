#ifdef USE_VTK
#include "set_vtk_data.h"
#include "global_value.h"
#include "reader_bin.h"

#include <filesystem>

namespace fs = std::filesystem;

int SetSolutionFromFileToVtk(const std::string &address_solution, vtkSmartPointer<vtkUnstructuredGrid> &u_grid, const bool sizeable) {

  std::vector<Type> data;
  std::vector<Type> size_data = {kRadiation, 1, 1, kDensity, kPressure};
  std::vector<std::string> name_data = {F_ILLUM, F_DIVSTREAM, F_ENERGY, F_DENSITY, F_PRESSURE};

  if (sizeable) {
    for (size_t i = 0; i < size_data.size(); i++)
      std::cout << "Size: " << name_data[i] << ": " << size_data[i] << std::endl;
  }

  std::vector<Vector3> data3;
  std::vector<std::string> name_data3 = {F_VELOCITY, F_STREAM, F_DIVIMPULS};

  std::vector<Matrix3> data9;
  std::vector<std::string> name_data9 = {F_IMPULS};

  const Type *sizes = size_data.data();
  for (auto &str : name_data) {
    if (fs::exists(address_solution + str)) {
      if (files_sys::bin::ReadSimple(address_solution + str, data)) {
        return e_completion_fail;
      }
      if (sizeable) {
        Type val = *sizes;
        for (size_t i = 0; i < data.size(); i++)
          data[i] *= val;
      }
      SetDoubleVtkData(fs::path(str).replace_extension(), data, u_grid);
    }
    sizes++;
  }

  for (auto &str : name_data3) {
    if (fs::exists(address_solution + str)) {
      if (files_sys::bin::ReadSimple(address_solution + str, data3)) {
        return e_completion_fail;
      }
      SetDoubleVtkData(fs::path(str).replace_extension(), data3, u_grid);
    }
  }

  for (auto &str : name_data9) {
    if (fs::exists(address_solution + str)) {
      if (files_sys::bin::ReadSimple(address_solution + str, data9)) {
        return e_completion_fail;
      }
      SetDoubleVtkData(fs::path(str).replace_extension(), data9, u_grid);
    }
  }

  return e_completion_success;
}
#endif //! USE_VTK