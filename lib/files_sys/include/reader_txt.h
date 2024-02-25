/**
 * @file reader_txt.h
 * @brief Чтение текстовых данных
 *
 */

#ifndef READER_TXT
#define READER_TXT

#include <fstream>
#include <set>
#include <string>
#include <vector>

#include "geo_types.h"
#include "solvers_struct.h"

#include "global_def.h"

/*! \addtogroup file_sys Файловый модуль
    @{
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
 * @return int ::e_type_completion
 */
template <typename Str_Type, typename T>
int ReadSimple(const Str_Type name_file, std::vector<T> &data) {
  std::ifstream ifile;
  OPEN_FSTREAM(ifile, std::string(name_file).c_str());

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
 * @brief Чтение текстового файла без размера данных
 *
 * @tparam Str_Type символьны тип
 * @tparam T тип считываемых данных
 * @param[in] name_file полное имя файла с расширением
 * @param[out] data набор std::set
 * @return int ::e_type_completion
 */
template <typename Str_Type, typename T>
int ReadSimple(const Str_Type name_file, std::set<T> &data) {
  std::ifstream ifile;
  OPEN_FSTREAM(ifile, std::string(name_file).c_str());

  int size;
  T buf;
  data.clear();

  ifile >> size;
  for (int i = 0; i < size; i++) {
    ifile >> buf;
    data.emplace(buf);
  }

  ifile.close();
  return e_completion_success;
}

/**
 * @brief Чтение текстового файла в набор уникальный ключей
 *
 * @tparam Str_Type Str_Type символьны тип
 * @tparam T тип считываемых данных
 * @param[in] name_file полное имя файла с расширением
 * @param[out] data массив std::vector
 * @return int ::e_type_completion
 */
template <typename Str_Type, typename T>
int ReadData(const Str_Type name_file, std::vector<T> &data) {

  std::ifstream ifile;
  OPEN_FSTREAM(ifile, name_file.c_str());
  int n = std::count(std::istreambuf_iterator<char>(ifile), std::istreambuf_iterator<char>(), '\n');
  ifile.close();

  data.resize(n);
  OPEN_FSTREAM(ifile, name_file.c_str());
  for (size_t i = 0; i < data.size(); i++) {
    ifile >> data[i];
  }

  ifile.close();
  return e_completion_success;
}

/**
 * @brief Чтение файла со сферой направлений
 *
 * @param[in] file_sphere_direction полное имя файла с расширением
 * @param[out] grid_direction структура сферы направлений
 * @note файл должен содержать размерность N, 3 координаты по каждому из N
 * направлений, дискретную площадь сферы
 * @return int ::e_type_completion
 */
int ReadSphereDirectionCartesian(const std::string &file_sphere_direction,
                                 grid_directions_t &grid_direction);

/**
 * @brief Чтение файла с координатами вершин и номерами граней внутренних ячеек
 *
 * @param[in] file_face_id полное имя файла с расширением
 * @param[out] inter_faces набор граней(std::map)
 * @warning ключом является номер ячейки, id - глобальная нумерация id_cell*cell_size + loc_id_face
 * @return int ::e_type_completion
 */
int ReadInitBoundarySetInFaces(const std::string &file_face_id, std::map<IntId, FaceCell> &inter_faces);

/**
 * @brief Чтение табулированной функции из файла
 *
 * @param file_func полное имя файла с расширением
 * @param tab_func структура табулированной функции
 * @return int ::e_type_completion
 */
int ReadTableFunc(const std::string &file_func, TableFunc &tab_func);
} // namespace txt
} // namespace files_sys

#endif // !READER_TXT
