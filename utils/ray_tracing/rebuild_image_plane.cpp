#include "ray_tracing_build_plane.h"
#include "utils.h"

int FUNC_NAME(RebuildImageBinToVtk)(int argc, char **argv) {
  if (argc != 7) {
    printf("Error input data!\n");
    printf("Input:\n");
    printf("path\\data.bin\n");
    printf("number_of_data\n");
    printf("plane_width\n");
    printf("plane_height\n");
    printf("number_of_pixels_width\n");
    printf("number_of_pixels_height\n");
    return e_completion_fail;
  }

  int size = std::stoi(argv[2]);
  Type X = std::stod(argv[3]);
  Type Y = std::stod(argv[4]);
  int X_pix = std::stoi(argv[5]);
  int Y_pix = std::stoi(argv[6]);
  ray_tracing::PlaneParams cfg(X,Y,X_pix,Y_pix);

  return ray_tracing::BuildVtkFromBin(cfg, size, argv[1]);
}
