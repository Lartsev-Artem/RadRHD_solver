#if 0 // def USE_VTK
#include "file_module.h"

#include "reader_vtk.h"
#include "writer_bin.h"

int FUNC_NAME(UnArchive)(int argc, char **argv) {
  utils::UtilsHelp(argc, argv, "Функция извлекает из архива файлы,"
                        "архив должен содержать файлы типа Solve0.vtk и Solve{last_idx}.vtk и файлы данных data.bin"
                        "затем вызывается процедура генерации сеток\n");

  std::string archive_name = ""; // "D:\\Desktop\\FilesCourse\\UtilsTest\\Solve.7z";

  switch (argc) {
  case 2:
    archive_name = argv[1];
    break;

  default:
    printf("Error input data!\n");
    printf("Input:\n");
    printf("path\\archive_name \n");
    printf("unzip directory\n\n");
    return 1;
  }

  std::filesystem::directory_entry fs_path(archive_name);
  std::string name = fs_path.path().filename().string();
  std::string path = std::regex_replace(archive_name, std::regex(name), "");

  std::string command = "C:\\\"Program Files\"\\7-Zip\\7z.exe x ";
  std::string cmd_call = (command + archive_name + " -o" + path);

  printf("zip cmd: %s\n", cmd_call.c_str());
  system(cmd_call.c_str());

  std::smatch pieces_match;
  std::regex_match(name, pieces_match, std::regex("(\\w+)\\.(\\w+)"));
  std::string new_path = path + pieces_match[1].str() + "\\"; //разархтвтрованный каталог

  std::string base_file;
  std::string base_name;

  int count_files = -1;
  for (const auto &entry : fs::directory_iterator(new_path)) //по всем файлам в директории
  {
    if (entry.is_regular_file()) {
      if (entry.path().extension() != ".vtk")
        continue;
      std::string file = entry.path().filename().string();
      std::regex_match(file, pieces_match, std::regex("(\\w+)\\.(\\w+)"));
      std::string clear_name = pieces_match[1].str();
      std::regex_search(clear_name, pieces_match, std::regex("\\s*(\\d+)"));

      if (std::stoi(pieces_match[0]) == 0) {
        base_file = entry.path().string();
        base_name = std::regex_replace(clear_name, std::regex("0"), "");
      }

      count_files = std::max(count_files, std::stoi(pieces_match[0]));
    }
  }

  char *new_argv[4];
  for (int i = 0; i < 4; i++)
    new_argv[i] = new char[200];

  // char new_argv[4][500];
  strcpy(new_argv[0], argv[0]);
  strcpy(new_argv[1], base_file.c_str());
  strcpy(new_argv[2], (new_path + base_name).c_str());
  strcpy(new_argv[3], std::to_string(count_files + 1).c_str());

  rebuild_solve(4, (char **)new_argv);

  for (int i = 0; i < 4; i++)
    delete[] new_argv[i];

  return e_completion_success;
}

#endif //! USE_VTK