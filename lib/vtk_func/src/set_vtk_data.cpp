#include "set_vtk_data.h"
#include "reader_bin.h"
#include <filesystem>
namespace fs = std::filesystem;

int SetSolutionFromFileToVtk(const std::string &address_solution, vtkSmartPointer<vtkUnstructuredGrid> &u_grid) {

  std::vector<Type> data;
  std::vector<std::string> name_data = {F_ILLUM, F_DIVSTREAM, F_ENERGY, F_DENSITY, F_PRESSURE};

  std::vector<Vector3> data3;
  std::vector<std::string> name_data3 = {F_VELOCITY, F_STREAM, F_DIVIMPULS};

  std::vector<Matrix3> data9;
  std::vector<std::string> name_data9 = {F_IMPULS};

  for (auto &str : name_data) {
    if (fs::exists(address_solution + str)) {
      if (files_sys::bin::ReadSimple(address_solution + str, data)) {
        return e_completion_fail;
      }
      SetDoubleVtkData(fs::path(str).replace_extension(), data, u_grid);
    }
  }

  for (auto &str : name_data3) {
    if (fs::exists(address_solution + str)) {
      if (files_sys::bin::ReadSimple(address_solution + str, data)) {
        return e_completion_fail;
      }
      SetDoubleVtkData(fs::path(str).replace_extension(), data, u_grid);
    }
  }

  for (auto &str : name_data9) {
    if (fs::exists(address_solution + str)) {
      if (files_sys::bin::ReadSimple(address_solution + str, data)) {
        return e_completion_fail;
      }
      SetDoubleVtkData(fs::path(str).replace_extension(), data, u_grid);
    }
  }

  return e_completion_success;
}