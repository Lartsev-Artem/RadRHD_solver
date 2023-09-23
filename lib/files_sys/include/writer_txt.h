/**
 * @file writer_txt.h
 * @brief Запись файлов в текстовом формате
 *
 */

#ifndef WRITER_TXT
#define WRITER_TXT

#include "dbgdef.h"
#include <fstream>
#include <typeinfo>
#include <vector>

/*! \addtogroup file_sys Файловый модуль
    @{
*/

namespace files_sys {
namespace txt {

/**
 * @brief Запись текстового файла
 *
 * @tparam Str_Type символьны тип
 * @tparam stdT std контейнер
 * @param[in] name_file полное имя файла с расширением
 * @param[in] data массив std::vector
 * @note  Файл содержит в первой строке число элементов. Далее последовательные
 * данные
 * @return size_t ::e_type_completion
 */
template <typename Str_Type, typename stdT>
size_t WriteSimple(const Str_Type name_file, const stdT &data) {
  std::ofstream ofile;
  OPEN_FSTREAM(ofile, std::string(name_file).c_str());

  ofile << data.size() << '\n';

  if (typeid(*data.begin()) == typeid(double)) { //если размер данных double, то выводить
    for (auto el : data)
      ofile << setprecision(16) << el << '\n';
  } else {
    for (auto el : data)
      ofile << el << '\n';
  }

  ofile.close();
  return e_completion_success;
}

} // namespace txt
} // namespace files_sys

#endif // !WRITER_TXT