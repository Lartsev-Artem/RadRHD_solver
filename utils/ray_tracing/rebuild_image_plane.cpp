#include "ray_tracing_build_plane.h"
#include "utils.h"

int FUNC_NAME(RebuildImageBinToVtk)(int argc, char **argv) {
  if (argc != 5) {
    printf("Error input data!\n");
    printf("Input:\n");
    printf("path\\data.bin\n");
    printf("number_of_data\n");
    printf("number_of_pixels_width\n");
    printf("number_of_pixels_height\n");
    return e_completion_fail;
  }

  int size = std::stoi(argv[2]);
  int X = std::stoi(argv[3]);
  int Y = std::stoi(argv[4]);
  return ray_tracing::BuildVtkFromBin(X, Y, size, argv[1]);
}
