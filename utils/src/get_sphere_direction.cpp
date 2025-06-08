#include "global_types.h"
#include "make_internal_format.h"
#include "utils.h"

int FUNC_NAME(GetSphereDirection)(int argc, char **argv)
{

  if (argc < 3)
  {
    printf("Error input data!\n");
    printf("Input: ");
    printf("path\\grid_sphere.vtk\n");
    printf("base_address\\ \n");
    printf("extended mode 0/1 [def=0]\\ \n");
    return e_completion_fail;
  }

  int mode = 0;
  if (argc == 4)
  {
    mode = std::stoi(argv[3]);

    if (mode != 0 && mode != 1)
    {
      printf("Error. Incorrect parameters mode=%d\n", mode);
      return e_completion_fail;
    }
  }

  return BuildSphereDirection(argv[1], argv[2], mode);
}