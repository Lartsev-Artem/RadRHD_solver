#include "file_module.h"

static int delete_files(const std::string &path, const std::string &mask) {
  const fs::path path_to_dir(path);
  const std::regex file_regex(mask);
  int count_del = 0;

  for (const auto &entry : fs::directory_iterator(path_to_dir)) //по всем файлам в директории
  {
    if (entry.is_regular_file()) {
      std::string file_name = entry.path().filename().string();

      if (std::regex_match(file_name, file_regex)) {
        printf("delete file: %s\n", entry.path().string().c_str());
        fs::remove(entry.path());
        count_del++;
      }
    }
  }

  return count_del;
}

int FUNC_NAME(DeleteSolveFiles)(int argc, char **argv) {
  if (argc != 5) {
    printf("Error input data!\n");
    printf("Input:\n");
    printf("path\n");
    printf("main_name_file (like \"Solve\")\n");
    printf("max number of files idx\n");
    printf("del mask\n");
    return 1;
  }

  const std::string path = argv[1];
  const std::string main_name = argv[2];
  const int n_files = std::stoi(argv[3]);
  const int del_mask = std::stoi(argv[4]);

  for (int i = 0; i < n_files; i++) //первый сохраняем
  {
    if (i % del_mask != 0) // оставить только каждый del_mask file (ex: 0,2,4,6,...)
    {
      std::string mask_file_regex(main_name + std::to_string(i) + "\[^0-9]*"); //все кроме цифр
      delete_files(path, mask_file_regex);
    }
  }

  return 0;
}
