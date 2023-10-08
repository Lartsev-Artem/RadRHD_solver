/**
 * @file file_module.h
 * @brief Функции оперирующие с файлами
 * @warning Могут быть не безопасны!!!

 */
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
/**
 * @brief Функция копирует группу файлов в указанную директория по маске
 *
 * @param argc[in] кол-во аргументов (требуется 5. 6 опционально)
 * @param argv[in] массив char** содержит
 * {src_path, dest_path, основное имя файла (обычно Solve),
 * максимальное число индексов файлов, доп.: шаг копирования}
 * @return int ::e_type_completion
 */
int CopySolveFiles(int argc, char **argv);

/**
 * @brief Функция удаляет группу файлов в указанную директория по маске
 *
 * @param argc[in] кол-во аргументов (требуется 5)
 * @param argv[in] массив char** содержит {
 * path, основное имя файла (обычно Solve),
 * максимальное число индексов файлов, шаг сохранения}
 * \note если шаг 4, то удалится: [save] 1,2,3 [save], 5,6,7 ...
 * @return int ::e_type_completion
 */
int DeleteSolveFiles(int argc, char **argv);

/**
 * @brief  Функция удаляет пропуски в индексации файлов
 *
 * @param argc[in] кол-во аргументов (требуется 4)
 * @param argv[in] массив char** содержит {path,
 * основное имя файла (обычно Solve),
 * максимальное число индексов файлов}
 * @return int ::e_type_completion
 */
int ReduceNameSolveFiles(int argc, char **argv);

/**
 * @brief Функция архивирует (вызовом 7z.exe) сетки vtk
 *
 * @details Функция разбирает набор сеток vtk на бинарные данные, пропуская первую
 * для сохранения геометрии и архивирует в единый архив. Опционально исходные сетки могут быть удалены
 * @param argc[in] кол-во аргументов (требуется 2. 3 опционально)
 * @param argv[in] массив char** содержит {path,основное имя файла (обычно Solve), признак удаления файлов}
 * @return int ::e_type_completion
 * @warning только под WINDOWS
 * @todo LINUX -tar
 */
int Archive(int argc, char **argv);

/**
 * @brief Функция из архива данных формирует набор сеток vtk
 *
 * @param argc[in] кол-во аргументов (требуется 2)
 * @param argv[in] массив char** содержит{имя архива, путь распаковки}
 * @return int ::e_type_completion
 */
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