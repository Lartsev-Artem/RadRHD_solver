
#include "bin_files_to_geo.h"
#include "utils.h"

int FUNC_NAME(ReWriteBinToGeo)(int argc, char **argv) {
  utils::UtilsHelp(argc, argv, "Функция на основе бинарных файлов геометрии формирует два файла с геометрией в формате ::face_t ::elem_t\n"
                               "на вход требуется указать путь к директории с геометрией");

  if (argc < 2) {
    RETURN_ERR("Need base address to geometry files\n");
  }

  return BinToGeo(std::string(argv[1]));
}