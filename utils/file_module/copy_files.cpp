#include "file_module.h"
#include "global_types.h"

static int copy_files(const std::string &src_path, const std::string &dest_path, const std::string &mask) {
  if (!fs::exists(dest_path)) {
    printf("no exist dest_dir: %s\n", dest_path.c_str());
    return 1;

    // fs::create_directory(dest_dir);
  }

  const fs::path path_to_dir(src_path);
  const std::regex file_regex(mask);
  int count_replace = 0;

  for (const auto &entry : fs::directory_iterator(path_to_dir)) //по всем файлам в директории
  {
    if (entry.is_regular_file()) {
      std::string file_name = entry.path().filename().string();

      if (std::regex_match(file_name, file_regex)) {
        fs::copy(entry.path(), dest_path, std::filesystem::copy_options::overwrite_existing);

        printf("copy file: %s\n", entry.path().string().c_str());
        count_replace++;
      }
    }
  }

  return count_replace;
}

int FUNC_NAME(CopySolveFiles)(int argc, char **argv) {
  if (argc < 5) {
    printf("Error input data!\n");
    printf("Input:\n");
    printf("src_path\n");
    printf("dest_path\n");
    printf("main_name_file (like \"Solve\")\n");
    printf("max number of files idx\n");
    printf("add: cpy mask\n");
    return e_completion_fail;
  }

  int cpy_mask = 1;
  if (argc == 6) {
    cpy_mask = std::stoi(argv[5]);
  }

  const std::string src_path = argv[1];
  const std::string dest_path = argv[2];
  const std::string main_name = argv[3];
  const int n_files = std::stoi(argv[4]);

  for (int i = 0; i < n_files; i += cpy_mask) {
    std::string mask_file_regex(main_name + std::to_string(i) + "\[^0-9]*");
    copy_files(src_path, dest_path, mask_file_regex);
  }

  return e_completion_success;
}
