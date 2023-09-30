#include "global_value.h"
#include "solvers_struct.h"

#include "graph_main.h"
#include "illum_main.h"
#include "trace_main.h"

#include "bin_files_to_geo.h"
#include "make_internal_format.h"
#include "trace_prebuild.h"

#include <filesystem>

namespace fs = std::filesystem;

int main(int argc, char **argv) {

  std::string prj_dir = "/home/artem/projects/solver/";

  glb_files.base_address = prj_dir + "tests/build/";
  glb_files.graph_address = prj_dir + "tests/build/graph/";
  glb_files.illum_geo_address = prj_dir + "tests/build/illum_geo/";
  glb_files.solve_address = prj_dir + "tests/build/Solve/";
  glb_files.name_file_sphere_direction = prj_dir + "tests/data/surface_26_dir.txt";
  glb_files.name_file_vtk = prj_dir + "tests/data/sphere.vtk";
  glb_files.Build();
  std::string file_vtk_out = glb_files.base_address + "illum_test.vtk";

  _solve_mode.class_vtk = e_grid_cfg_default;
  _solve_mode.accuracy = 1e-5;
  _solve_mode.max_number_of_iter = 1;

  fs::create_directory(glb_files.base_address);
  fs::create_directory(glb_files.graph_address);
  fs::create_directory(glb_files.illum_geo_address);
  fs::create_directory(glb_files.solve_address);
  fs::copy_file(glb_files.name_file_vtk, file_vtk_out, fs::copy_options::overwrite_existing);
  glb_files.solve_address += "Solve";

  if (BuildDataFromVTK(glb_files) != e_completion_success) {
    RETURN_ERR("don't build start format\n");
  }

  if (trace::PreBuild(glb_files) != e_completion_success) {
    RETURN_ERR("don't build trace struct\n");
  }

  if (BinToGeo(glb_files.base_address) != e_completion_success) {
    RETURN_ERR("don't build geo format\n");
  }

  graph::RunGraphModule();
  trace::RunTracesModule();
  illum::RunIllumModule();

  return 0;
}