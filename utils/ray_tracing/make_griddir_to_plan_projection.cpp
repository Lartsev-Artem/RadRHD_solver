#include "utils.h"

#include "writer_txt.h"
#include "writer_bin.h"
#include "global_value.h"

#include "ray_tracing_build_plane.h"

int FUNC_NAME(MakeGridDirToPlanProjection)(int argc, char **argv) {
    if (argc != 12) {
    printf("Error input data!\n");
    printf("Input:\n");
    printf("path\\output_dir\n");
    printf("observer's position: X Y Z\n");
    printf("plane origin: X Y Z\n");
    printf("plane width\n");
    printf("plane height\n");
    printf("number_of_pixels_width\n");
    printf("number_of_pixels_height\n");
    return e_completion_fail;
  }
  
  std::string output_dir = argv[1];
  Vector3 observer(std::stod(argv[2]), std::stod(argv[3]), std::stod(argv[4]));

  Vector3 plane_orig(std::stod(argv[5]), std::stod(argv[6]), std::stod(argv[7]));
  Type plane_width = std::stod(argv[8]);
  Type plane_height = std::stod(argv[9]);

  int pixels_width = std::stoi(argv[10]);
  int pixels_height = std::stoi(argv[11]);
  
  ray_tracing::PlaneParams plane(plane_width,plane_height,pixels_width,pixels_height);
  ray_tracing::ParamTraceProjection params(plane,plane_orig,observer);

  std::vector<Ray_t> rays;    
  ray_tracing::MakeRays(plane,plane_orig,observer,rays);

  grid_directions_t grid_dir(rays.size());

    for (size_t i = 0; i < grid_dir.size; i++)
    {        
        grid_dir.directions[i].area = 1;
        grid_dir.directions[i].dir = -rays[i].direction; //у единичной сферы начало направления- точка (0,0,0)        
    }
    grid_dir.full_area = grid_dir.size;


  if (files_sys::bin::WriteSimple (output_dir + F_PLANE_CFG, sizeof(params) , (uint8_t*)(&params)))
    return e_completion_fail;

  if (files_sys::txt::WriteSphereDirectionCartesian(output_dir+F_DIRECTION_GRID, grid_dir))
    return e_completion_fail;

  return e_completion_success;
}
