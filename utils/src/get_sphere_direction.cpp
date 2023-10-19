#include "global_types.h"
#include "make_internal_format.h"
#include "utils.h"

int FUNC_NAME(GetSphereDirection)(int argc, char **argv) {

  if (argc != 3) {
    printf("Error input data!\n");
    printf("Input: ");
    printf("path\\grid_sphere.vtk\n");
    printf("base_address\\ \n");
    return e_completion_fail;
  }

  return BuildSphereDirection(argv[1], argv[2]);
}