#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <string>

//#define GENERATE_UTILS  ///<генерировать исполняемые файлы из утилит

#ifdef GENERATE_UTILS
#define FUNC_NAME(_name) main
#else
#define FUNC_NAME(_name) utils::_name
#endif

namespace utils {

template <typename T, typename T2>
void UtilsHelp(int argc, T str, T2 msg = "") {
  if (argc == 1) {
    return; //нет аргументов
  }

  std::string word = str[1];

  if (word.find("-h", 0, 2) == std::string::npos && word.find("--help", 0) == std::string::npos) {
    return;
  }

  setlocale(LC_CTYPE, "rus");
  std::cout << msg << '\n';
  exit(0);
}

#ifndef GENERATE_UTILS
int ReBuildNetgenToVTK(int argc, char **argv);
int SetScalarDataVtkFromBinFile(int argc, char **argv);
int DataToSimpleTxt(int argc, char **argv);
int AddScalarDataVtkFromBinFile(int argc, char **argv);
#endif

} // namespace utils

#endif //! UTILS_H
