/**
 * @file writer_txt.h
 * @brief Запись файлов в текстовом формате
 *
 */

#ifndef WRITER_TXT
#define WRITER_TXT

#include <fstream>
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
 * @tparam T тип записываемых данных
 * @param[in] name_file полное имя файла с расширением
 * @param[in] data массив std::vector
 * @note  Файл содержит в первой строке число элементов. Далее последовательные
 * данные
 * @return size_t ::e_type_completion
 */
template <typename Str_Type, typename T>
size_t WriteSimple(const Str_Type name_file, std::vector<T> &data) {
  std::ofstream ofile;
  OPEN_FSTREAM(ofile, name_file.c_str());

  ofile << data.size() << '\n';

  for (int i = 0; i < data.size(); i++) {
    ofile << data[i] << '\n';
  }

  ofile.close();
  return e_completion_success;
}

} // namespace txt
} // namespace files_sys

#endif // !WRITER_TXT