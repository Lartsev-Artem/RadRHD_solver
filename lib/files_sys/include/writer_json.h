/**
 * @file writer_json.h
 * @brief Файл содержит шаблон для записи json структуры в файл
 *
 */

#ifndef WRITER_JSON
#define WRITER_JSON

#include <fstream>
#include <string>

#include "dbgdef.h"
#include "global_types.h"
#include "json/json.hpp"


/*! \addtogroup file_sys Файловый модуль
    @{
*/

namespace files_sys {
namespace json {

/**
 * @brief
 *
 * @tparam file_t строка
 * @tparam struct_t тип json структуры
 * @param[in] file полное имя файла с расширением
 * @param[in] st json структура
 * @param[in] mod форма записи
 * @return int ::e_type_completion
 */
template <typename file_t, typename struct_t>
int WriteFileJson(const file_t file, struct_t &st, std::ios_base::openmode mod = std::ios::trunc) {
  std::ofstream ofile(file, mod);

  if (!ofile.is_open()) {
    RETURN_ERR("Error : file %s is not open\n", std::string(file).c_str());
  }

  nlohmann::json j;
  j = st;

  ofile.width(1);
  ofile << j;
  ofile.close();
  return e_completion_success;
}

} // namespace json
} // namespace files_sys
#endif //! WRITER_JSON