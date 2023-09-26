#ifdef USE_VTK
#include "file_module.h"

#include "reader_vtk.h"
#include "writer_bin.h"

int FUNC_NAME(Archive)(int argc, char **argv) {
  utils::UtilsHelp(argc, argv, "Функция извлекает из последовательности сеток файлы с решением (Solvei.vtk),"
                               "удаляет существующие файлы сеток(если установлен флаг delete)\n"
                               "сохраняет начальную и финальную конфигурацию");

  std::string files_path; // = "D:\\Desktop\\FilesCourse\\UtilsTest\\";
  std::string files_name; // = "Solve";
  bool is_delete = false;

  switch (argc) {
  case 3:
    files_path = argv[1];
    files_name = argv[2];
    break;
  case 4:
    files_path = argv[1];
    files_name = argv[2];
    is_delete = std::stoi(argv[3]);
    break;

  default:
    printf("Error input data!\n");
    printf("Input:\n");
    printf("full_path \n");
    printf("NameGrid (e.g. Solve, if files (Solve0.vtk,...,SolveN.vtk))\n");
    printf("add param: delete (1 - delete base files ), 0 - don't do this (default)\n");
    return e_completion_fail;
  }

  if (!fs::exists(files_path)) {
    printf("no exist files_dir: %s\n", files_path.c_str());
    return e_completion_fail;
  }

  const std::regex file_regex(files_name);
  int count_files = 0;
  std::string file_name = "";

  for (const auto &entry : fs::directory_iterator(files_path)) //по всем файлам в директории
  {
    if (entry.is_regular_file()) {
      if (entry.path().extension() != ".vtk")
        continue;

      if (count_files > 1 && is_delete) {
        std::remove((files_path + file_name).c_str());
      }
      file_name = entry.path().filename().string();
      if (std::regex_search(file_name, file_regex)) {
        vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

        if (files_sys::vtk::Read(entry.path().string(), ugrid))
          RETURN_ERR("Error reading the file vtk\n");

        const int size_grid = ugrid->GetNumberOfCells();

        if (ugrid->GetCellData()->GetNumberOfArrays() == 0)
          RETURN_ERR("grid hasn't data\n");

        std::vector<Type> v_data;
        std::vector<Vector3> v_data3;
        std::vector<Matrix3> v_data9;
        Type *v;

        for (size_t i = 0; i < ugrid->GetCellData()->GetNumberOfArrays(); i++) {
          std::string name_data(ugrid->GetCellData()->GetArrayName(i));
          vtkDataArray *data = ugrid->GetCellData()->GetScalars(name_data.c_str());

          int N = data->GetNumberOfComponents();
          std::string file_bin = files_path + files_name + std::to_string(count_files) + name_data + ".bin";

          switch (N) {
          case 1:
            v_data.resize(size_grid);
            for (size_t i = 0; i < size_grid; i++) {
              v_data[i] = data->GetTuple1(i);
            }
            files_sys::bin::WriteSimple(file_bin, v_data);
            break;

          case 3:
            data = ugrid->GetCellData()->GetVectors(name_data.c_str());
            v_data3.resize(size_grid);
            for (size_t i = 0; i < size_grid; i++) {
              v = data->GetTuple3(i);
              for (size_t j = 0; j < 3; j++) {
                v_data3[i](j) = v[j];
              }
            }
            files_sys::bin::WriteSimple(file_bin, v_data3);
            break;

          case 9:
            data = ugrid->GetCellData()->GetTensors(name_data.c_str());
            v_data9.resize(size_grid);
            for (size_t i = 0; i < size_grid; i++) {
              v = data->GetTuple9(i);
              for (size_t j = 0; j < 9; j++) {
                v_data9[i](j) = v[j];
              }
            }
            files_sys::bin::WriteSimple(file_bin, v_data9);
            break;
          default:
            RETURN_ERR((std::string("unknow data type") + name_data).c_str());
          }
        }

        count_files++;
      }
    }
  }

  std::string command = "C:\\\"Program Files\"\\7-Zip\\7z.exe a";
  std::string zip_key = " -mx9 -sdel "; //степень сжатия, удалять файлы после создания
  std::string name_archive = files_path + files_name + ".7z ";
  std::string cmd_call = (command + zip_key + name_archive + files_path);

  printf("zip cmd: %s\n", cmd_call.c_str());
  system(cmd_call.c_str());

  return e_completion_success;
}
#endif //! USE_VTK