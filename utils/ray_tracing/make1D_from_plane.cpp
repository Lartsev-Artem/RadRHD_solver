#include "ray_tracing_build_plane.h"
#include "utils.h"

#include "reader_vtk.h"
#include "writer_txt.h"

int FUNC_NAME(Make1DFromPlane)(int argc, char **argv) {
  if (argc != 11) {
    printf("Error input data!\n");
    printf("Input:\n");
    printf("path\\file_grid.vtk\n");
    printf("path\\output_trace.txt\n");
    printf("plane_size: width, height\n");    
    printf("number_of_pixels: width, height\n");    
    printf("start ray: X, Y\n");
    printf("direction: X, Y\n");    
    return e_completion_fail;
  }
  
  std::string file_grid = argv[1];
  std::string file_trace = argv[2];

  Type X = std::stod(argv[3]);
  Type Y = std::stod(argv[4]);
  int X_pix = std::stoi(argv[5]);
  int Y_pix = std::stoi(argv[6]);
      
  Vector2 orig = Vector2(std::stod(argv[7]),std::stod(argv[8]));
  Eigen::Vector2i dir = Eigen::Vector2i(std::stoi(argv[9]),std::stoi(argv[10]));

  Vector3 angle_of_plane = Vector3(-(X / 2), -(Y / 2), 0);
  Type step_x = X / X_pix;
  Type step_y = Y / Y_pix;
  
  //позиция центрального пикселя
  int i_px = (orig[0]-angle_of_plane[0]) / step_x; 
  int i_py = (orig[1]-angle_of_plane[1]) / step_y;

  std::vector<Type> data;
  vtkSmartPointer<vtkUnstructuredGrid> grid;
  if(files_sys::vtk::Read(file_grid,grid))
  {
    printf("File %s not opened!\n", file_grid.c_str());
    return e_completion_fail;
  }

  if(grid->GetNumberOfCells() != X_pix*Y_pix)
  {
    printf("Error size plane: grid:%lld, Plane: %dx%d\n", grid->GetNumberOfCells(), X_pix,Y_pix);
    return e_completion_fail;
  }

  if(dir == Eigen::Vector2i(1,0))
  {
    data.resize(X_pix-i_px);
    for (size_t i = i_px; i < X_pix; i++)
    {      
      data[i-i_px] = grid->GetCellData()->GetArray(0)->GetTuple1(i*Y_pix + i_py);      
    }    
  }
  else if(dir == Eigen::Vector2i(-1,0))
  {    
     data.resize(i_px-0);
    for (size_t i = i_px; i >= 0; i--)
    {      
      data[i] = grid->GetCellData()->GetArray(0)->GetTuple1(i*Y_pix + i_py);      
    } 

  }
  else if(dir == Eigen::Vector2i(0,1))
  {
     data.resize(Y_pix-i_py);
    for (size_t i = i_py; i < Y_pix; i++)
    {      
      data[i-i_py] = grid->GetCellData()->GetArray(0)->GetTuple1(i_px*Y_pix + i);      
    }  
  }
  else if(dir == Eigen::Vector2i(0,-1))
  {
    data.resize(i_py - 0);
    for (size_t i = i_py; i >= 0; i--)
    {      
      data[i] = grid->GetCellData()->GetArray(0)->GetTuple1(i_px*Y_pix + i);      
    }      
  }
  else if(dir == Eigen::Vector2i(1,1))
  {
    data.resize(X_pix*Y_pix);
    std::vector<Type> R(X_pix*Y_pix);
    ray_tracing::PlaneParams cfg(X,Y,X_pix,Y_pix);
    for (size_t i = 0; i < X_pix; i++)
    for (size_t j = 0; j < Y_pix; j++)
    {      
      R[i*Y_pix + j] = cfg.get_pixel_coord(i,j).norm();
      data[i*Y_pix + j] = grid->GetCellData()->GetArray(0)->GetTuple1(i*Y_pix + j);
    }  

    files_sys::txt::WriteSimple(file_trace,R, data);
    return e_completion_success;
  }
  else
  {
    printf("Error: Unavailable direction. Use only (+-1,0), (0, +-1) or (1,1) to radius projection\n");
    return e_completion_fail;
  }

  files_sys::txt::WriteSimple(file_trace, data);
  

  return e_completion_success;
}
