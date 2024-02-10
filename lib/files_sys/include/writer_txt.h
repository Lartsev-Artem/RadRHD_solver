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

#include "geo_types.h"

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
 * @return int ::e_type_completion
 */
template <typename Str_Type, typename stdT>
int WriteSimple(const Str_Type name_file, const stdT &data) {
  std::ofstream ofile;
  OPEN_FSTREAM(ofile, std::string(name_file).c_str());

  ofile << data.size() << '\n';

  if (typeid(*data.begin()) == typeid(double)) { //если размер данных double, то выводить
    for (auto el : data)
      ofile << std::setprecision(16) << el << '\n';
  } else {
    for (auto el : data)
      ofile << el << '\n';
  }

  ofile.close();
  return e_completion_success;
}

template <typename Str_Type, typename stdT>
int WriteSimple(const Str_Type name_file, const stdT &x, const stdT &y) {
  std::ofstream ofile;
  OPEN_FSTREAM(ofile, std::string(name_file).c_str());

  if (x.size() != y.size()) {
    RETURN_ERR("size x!= size y\n");
  }

  ofile << x.size() << '\n';

  auto it = y.data();
  if (typeid(*x.begin()) == typeid(double)) { //если размер данных double, то выводить
    for (auto el : x) {
      ofile << std::setprecision(16) << el << ' ' << *it << '\n';
      it++;
    }
  } else {
    for (auto el : y)
      ofile << el << ' ' << *it << '\n';
  }

  ofile.close();
  return e_completion_success;
}

/**
 * @brief
 *
 * @param[in] file_sphere_direction полное имя файла с расширением
 * @param[in] grid_direction сетка сферы направлений
 * @return int ::e_type_completion
 */
int WriteSphereDirectionCartesian(const std::string &file_sphere_direction, const grid_directions_t &grid_direction);

} // namespace txt
} // namespace files_sys

#endif // !WRITER_TXT