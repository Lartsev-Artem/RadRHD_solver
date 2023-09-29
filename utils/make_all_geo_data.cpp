#ifdef USE_VTK
#include "utils.h"

#include <filesystem>

#include "bin_files_to_geo.h"
#include "make_internal_format.h"
#include "trace_prebuild.h"

namespace fs = std::filesystem;

int FUNC_NAME(MakeAllGeoData)(int argc, char **argv) {
  utils::UtilsHelp(argc, argv, "Функция строит для всех модулей геометрические настроечные файлы из vtk сетки\n");

  if (argc < 3) {
    RETURN_ERR("Need base_address and file_vtk\n");
  }

  global_files_t files;
  files.base_address = argv[1];
  files.name_file_vtk = argv[2];
  files.Build();

  if (!fs::exists(files.base_address) || !fs::exists(files.name_file_vtk)) {
    RETURN_ERR("bad input address\n");
  }

  if (BuildDataFromVTK(files) != e_completion_success) {
    return e_completion_fail;
  }

  if (trace::PreBuild(files) != e_completion_success) {
    return e_completion_fail;
  }

  if (BinToGeo(files.base_address) != e_completion_success) {
    return e_completion_fail;
  }

  return e_completion_success;
}
#endif