#include "utils.h"

#include "ray_tracing_build_plane.h"

#include "writer_vtk.h"
#include "reader_bin.h"

#include <vtkDoubleArray.h>

int FUNC_NAME(MakeObserverPlane)(int argc, char **argv) {
    if (argc != 5) {
    printf("Error input data!\n");
    printf("Input:\n");
    printf("path\\Data_i.bin\n");
    printf("path\\projection_file.bin\n");
    printf("path\\plane_cfg_file.bin\n");
    printf("path\\image_plane.vtk\n");    
    return e_completion_fail;
  }
  
  std::string address_solve = argv[1];      
  std::string file_projection = argv[2];
  std::string file_plane_cfg = argv[3];
  std::string output_grid = argv[4];

  ray_tracing::ParamTraceProjection plane_cfg(ray_tracing::PlaneParams(1,1,1,1),Vector3::Zero(),Vector3::Zero());
  if (files_sys::bin::ReadSimple(file_plane_cfg, (uint8_t*)(&plane_cfg)))
    return e_completion_fail;

  ray_tracing::PlaneParams plane = plane_cfg.params2D;

  vtkSmartPointer<vtkUnstructuredGrid> image_plane = vtkSmartPointer<vtkUnstructuredGrid>::New();
  
  ray_tracing::MakeVtkPlane(plane,image_plane);
  
  // проверка наличия файлов данных

  std::vector<IntId> projection;
  if (files_sys::bin::ReadSimple(file_projection, projection)) {
    return e_completion_fail;
  }

  vtkSmartPointer<vtkDoubleArray> data_vtk = vtkSmartPointer<vtkDoubleArray>::New();
  data_vtk->SetNumberOfComponents(1);
  data_vtk->SetName("data");
  std::vector<Type> data;
  for (size_t i = 0; i < image_plane->GetNumberOfCells(); i++)
  {    
     if(files_sys::bin::ReadSimple(address_solve+std::to_string(i)+".bin",data))
     {
        printf("Error reading\n");
        return e_completion_fail;
     }

      int id_cell = projection[i];
      if (id_cell < (int)data.size() * CELL_SIZE) 
      {
        double val = id_cell < 0 ? 0 : data[projection[i]/CELL_SIZE];
        data_vtk->InsertNextTuple(&val);        
      } 
      else 
      {
        std::cout << "err projection. Id= " << id_cell << ", N=" << data.size() << "\n";
        return e_completion_fail;
      }     
  }
  image_plane->GetCellData()->AddArray(data_vtk);

  if (files_sys::vtk::WriteVtkGrid(output_grid, image_plane))
      return e_completion_fail;

  return e_completion_success;
}
