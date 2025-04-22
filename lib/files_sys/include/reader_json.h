/**
 * @file reader_json.h
 * @brief Чтение json данных
 *
 */

#ifndef READER_JSON
#define READER_JSON

#include <string>

#include "global_value.h"
#include "json_struct.h"

/*! \addtogroup file_sys Файловый модуль
    @{
*/
namespace files_sys {

/**
 * @brief Пространство имён подмодуля json файлов
 *
 */
namespace json {

/**
 * @brief Чтение файла со структурой проекта (путей  файлов)
 *
 * @param[in] file полное имя файла с расширением
 * @param[out] st структура путей
 * @return int ::e_type_completion
 */
int Read(const std::string &file, global_files_t &st);

/**
 * @brief Чтение файла настройкой глобального решателя
 *
 * @param[in] file полное имя файла с расширением
 * @param[out] st структура данных
 * @return int ::e_type_completion
 */
int Read(const std::string &file, solve_mode_t &st);

/**
 * @brief Чтение файла настройкой газодинамического решателя
 *
 * @param[in] file полное имя файла с расширением
 * @param[out] st структура данных
 * @return int ::e_type_completion
 */
int Read(const std::string &file, hllc_value_t &st);

/**
 * @brief Чтение общей конфигурации проекта
 *
 * @tparam file символный тип данных
 * @param[in] file_set полное имя файла с расширением
 * @param[out] glb_files структура путей
 * @param[out] solve_mode структура настройки глобального решателя
 * @param[out] hllc_conf структура настройки газодинамического решателя
 * @return int ::e_type_completion
 */
template <typename file>
int ReadStartSettings(file file_set, global_files_t &glb_files,
                      solve_mode_t *solve_mode = nullptr,
                      hllc_value_t *hllc_conf = nullptr) {
  if (Read(file_set, glb_files)) {
    return e_completion_fail;
  }

  glb_files.name_file_settings = file_set;

  glb_files.Build();

#ifdef DEBUG
  glb_files.print();
#endif

  if (solve_mode != nullptr) {
    if (Read(glb_files.solve_configuration, *solve_mode))
      return e_completion_fail;
  }

  if (hllc_conf != nullptr) {
    if (Read(glb_files.name_file_hllc_set, *hllc_conf))
      return e_completion_fail;
  }
  return e_completion_success;
}

} // namespace json
} // namespace files_sys

#endif // !READER_JSON
