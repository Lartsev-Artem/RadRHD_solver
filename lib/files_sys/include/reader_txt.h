/**
 * @file reader_txt.h
 * @brief Чтение текстовых данных
 *
 */

#ifndef READER_TXT
#define READER_TXT

#include <fstream>
#include <string>
#include <vector>

#include "../global_def.h"
#include "geo_types.h"

/*! \addtogroup file_sys Файловый модуль
    @{
*/

/**
 * @brief Пространство имён файлового модуля
 *
 */
namespace files_sys {

/**
 * @brief Пространство имён подмодуля текстовых файлов
 *
 */
namespace txt {

/**
 * @brief Чтение текстового файла
 *
 * @tparam Str_Type символьны тип
 * @tparam T тип считываемых данных
 * @param[in] name_file полное имя файла с расширением
 * @param[out] data массив std::vector
 * @note Файл должен содержать в первой строке число элементов. Далее
 * последовательные данные.
 * @return size_t ::e_type_completion
 */
template <typename Str_Type, typename T>
size_t ReadSimple(const Str_Type name_file, std::vector<T> &data) {
  std::ifstream ifile;
  OPEN_FSTREAM(ifile, name_file);

  int size;
  ifile >> size;
  data.resize(size);

  for (int i = 0; i < size; i++) {
    ifile >> data[i];
  }

  ifile.close();
  return e_completion_success;
}

/**
 * @brief Чтение файла со сферой направлений
 *
 * @param[in] file_sphere_direction
 * @param[out] grid_direction
 * @note файл должен содержать размерность N, 3 координаты по каждому из N
 * направлений, дискретную площадь сферы
 * @return int ::e_type_completion
 */
int ReadSphereDirectionСartesian(const std::string &file_sphere_direction,
                                 grid_directions_t &grid_direction);

} // namespace txt
} // namespace files_sys

#endif // !READER_TXT
