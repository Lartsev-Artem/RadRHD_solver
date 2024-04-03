#include "ray_tracing_build_plane.h"
#include "utils.h"

#include "reader_bin.h"
#include "reader_vtk.h"
#include "writer_vtk.h"

#include "set_vtk_data.h"
#include <vtkDoubleArray.h>

static int get_data_projection(const std::vector<IntId> &projection, vtkDataArray *data, std::vector<Type> &s_data) {

  int components = data->GetNumberOfComponents();
  int size = (data->GetSize() - components) / components;
  switch (components) {
  case 1: {
    for (size_t p = 0; p < s_data.size(); p++) {
      int id_cell = projection[p];
      if (id_cell < size * CELL_SIZE) {
        s_data[p] = id_cell < 0 ? 0 : data->GetTuple1(id_cell / CELL_SIZE);
      } else {
        std::cout << "err projection. Id= " << id_cell << ", N=" << size << "\n";
        return e_completion_fail;
      }
    }
    break;
  }
  case 3: {

    Vector3 buf;
    for (size_t p = 0; p < s_data.size(); p++) {
      int id_cell = projection[p];
      if (id_cell < size * CELL_SIZE) {
        for (int j = 0; j < components; j++) {
          buf[j] = id_cell < 0 ? 0 : data->GetTuple3(id_cell / CELL_SIZE)[j];
        }
        s_data[p] = buf.norm();
      } else {
        std::cout << "err projection. Id= " << id_cell << ", N=" << size << "\n";
        return e_completion_fail;
      }
    }
    break;
  }
  case 9: {
    Matrix3 buf;
    for (size_t p = 0; p < s_data.size(); p++) {
      int id_cell = projection[p];
      if (id_cell < size * CELL_SIZE) {
        for (int j = 0; j < components; j++) {
          buf.data()[j] = id_cell < 0 ? 0 : data->GetTuple9(id_cell / CELL_SIZE)[j];
        }
        s_data[p] = buf.norm();
      } else {
        std::cout << "err projection. Id= " << id_cell << ", N=" << size << "\n";
        return e_completion_fail;
      }
    }
    break;
  }

  default:
    s_data.assign(s_data.size(), 0);
    break;
  }

  return e_completion_success;
}

int FUNC_NAME(PlaneFromVtk)(int argc, char **argv) {
  if (argc != 7) {
    printf("Error input data!\n");
    printf("Input:\n");
    printf("path\\Solve\n");
    printf("path\\projection_id_cell.bin\n");
    printf("number_of_data\n");
    printf("number_of_pixels_width\n");
    printf("number_of_pixels_height\n");
    printf("path\\Output\n");
    return e_completion_fail;
  }

  const std::string address_solve = argv[1];
  const std::string file_projection = argv[2];
  int number_of_data = std::stoi(argv[3]);
  int X = std::stoi(argv[4]);
  int Y = std::stoi(argv[5]);
  const std::string output = argv[6];

  std::vector<IntId> projection;
  if (files_sys::bin::ReadSimple(file_projection, projection)) {
    return e_completion_fail;
  }

  if (X * Y != projection.size()) {
    std::cout << "X*Y != projection \n";
    return e_completion_fail;
  }

  std::vector<Type> s_data(projection.size());
  vtkSmartPointer<vtkDoubleArray> data_array = vtkSmartPointer<vtkDoubleArray>::New();

  for (int num = 0; num < number_of_data; num++) {

    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    if (files_sys::vtk::Read(address_solve + std::to_string(num) + ".vtk", grid)) {
      return e_completion_fail;
    }

    // make plane
    vtkSmartPointer<vtkUnstructuredGrid> grid_plane = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ray_tracing::MakeVtkPlane(X, Y, grid_plane);

    for (size_t i = 0; i < grid->GetCellData()->GetNumberOfArrays(); i++) {
      std::string name_data(grid->GetCellData()->GetArrayName(i));
      vtkDataArray *data = grid->GetCellData()->GetScalars(name_data.c_str());
      if (get_data_projection(projection, data, s_data) == e_completion_success) {
        // add data to plane
        SetDoubleVtkData(name_data, s_data, grid_plane);
      }
    }

    if (files_sys::vtk::WriteVtkGrid(output + std::to_string(num) + ".vtk", grid_plane))
      return e_completion_fail;
  }

  return e_completion_success;
}
