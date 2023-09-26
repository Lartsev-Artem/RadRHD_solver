#ifndef FILE_MODULE_H
#define FILE_MODULE_H

#include "utils.h"

#include <fstream>
#include <iostream>
#include <string>

#include <filesystem> // need -std::c++17
#include <regex>      //регулярные выражения

#if _MSC_VER
#define VERSION _MSVC_LANG
#else
#define VERSION __cplusplus
#endif

#if VERSION < 201703
#error "Need c++17 to use file_utils"
#endif

namespace fs = std::filesystem;

#ifndef GENERATE_UTILS
namespace utils {
int CopySolveFiles(int argc, char **argv);
int DeleteSolveFiles(int argc, char **argv);
int ReduceNameSolveFiles(int argc, char **argv);
int Archive(int argc, char **argv);
int UnArchive(int argc, char **argv);
} // namespace utils
#endif

#if VERSION < 201703
inline char *_strdup(const char *s) {
  char *val = (char *)malloc(sizeof(s));
  return strcpy(val, s);
}
#else
#define _strdup(s) strdup(s);
#endif

#endif //! FILE_MODULE_H