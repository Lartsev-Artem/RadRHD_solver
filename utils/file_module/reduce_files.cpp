#include "file_module.h"

static bool rename_files(const std::string &path, const std::string &mask, const int cur_idx, const int new_idx) {
  const fs::path path_to_dir(path);
  const std::regex file_regex(mask);
  const std::regex replace_regex(std::to_string(cur_idx));

  bool exist_f = false;

  for (const auto &entry : fs::directory_iterator(path_to_dir)) //по всем файлам в директории
  {
    if (entry.is_regular_file()) {
      std::string file_name = entry.path().filename().string();

      if (std::regex_match(file_name, file_regex)) {
        std::string path = entry.path().generic_string();
        path.erase(path.end() - file_name.length(), path.end()); // путь к файлу без имени

        std::string new_file_name = std::regex_replace(file_name, replace_regex, std::to_string(new_idx));
        fs::rename(entry.path(), path + new_file_name);

        exist_f = true;

        printf("rename: %s -> %s\n", file_name.c_str(), new_file_name.c_str());
      }
    }
  }

  return exist_f;
}

int FUNC_NAME(ReduceNameSolveFiles)(int argc, char **argv) {
  if (argc != 4) {
    printf("Error input data!\n");
    printf("Input:\n");
    printf("path\n");
    printf("main_name_file (like \"Solve\")\n");
    printf("max number of files idx\n");
    return 1;
  }

  const std::string path = argv[1];
  const std::string main_name = argv[2];
  const int n_files = std::stoi(argv[3]);

  int new_idx = 0;
  for (int i = 0; i < n_files; i++) {
    std::string mask_file_regex(main_name + std::to_string(i) + "\[^0-9]*");

    if (rename_files(path, mask_file_regex, i, new_idx)) {
      new_idx++;
    }
  }
}
