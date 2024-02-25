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

  fs::create_directory(files.base_address + "trace");
  fs::create_directory(files.base_address + "graph");
  fs::create_directory(files.base_address + "illum_geo");
  fs::create_directory(files.base_address + "Solve");

  if (!fs::exists(files.base_address) || !fs::exists(files.name_file_vtk)) {
    RETURN_ERR("bad input address\n");
  }

  if (BuildDataFromVTK(files) != e_completion_success) {
    RETURN_ERR("don't build start format\n");
  }

  if (trace::PreBuild(files) != e_completion_success) {
    RETURN_ERR("don't build trace struct\n");
  }

  if (BinToGeo(files.base_address) != e_completion_success) {
    RETURN_ERR("don't build geo format\n");
  }

  fs::copy_file(files.name_file_vtk, files.base_address);

  return e_completion_success;
}
#endif