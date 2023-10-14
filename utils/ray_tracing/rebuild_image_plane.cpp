#include "ray_tracing_build_plane.h"
#include "utils.h"

int FUNC_NAME(RebuildImageBinToVtk)(int argc, char **argv) {
  if (argc != 3) {
    printf("Error input data!\n");
    printf("Input:\n");
    printf("path\\data.bin\n");
    printf("number_of_data\n");
    return e_completion_fail;
  }

  int size = std::stoi(argv[2]);
  return ray_tracing::BuildVtkFromBin(size, argv[1]);
}
