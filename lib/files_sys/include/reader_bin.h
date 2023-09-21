#ifndef READER_BIN
#define READER_BIN

#include <cstdio>
#include <vector>

#include "dbgdef.h"
#include "geo_types.h"

namespace files_sys {
namespace bin {

/**
 * @brief Чтение бинарного файла
 *
 * @tparam T - тип данных
 * @param[in] name_file полное имя файла с расширением
 * @param[out] data массив std::vector
 * @note Файл должен содержать в первой строке число элементов. Далее
 * последовательные данные.
 * @warning размер типа данных должен соответствовать данным в файле
 * @return size_t - ::e_type_completion
 */
template <typename T>
size_t ReadSimple(const std::string& name_file, std::vector<T>& data) {
  FILE* f;
  OPEN_FILE(f, name_file.c_str(), "rb");

  int n;
  fread(&n, sizeof(int), 1, f);
  data.resize(n);
  fread(data.data(), sizeof(T), n, f);

  fclose(f);
  return e_completion_success;
}

/**
 * @brief Чтение геометрии элементов сетки
 *
 * @tparam geo_elem - структура, содержащая поле "geo"
 * @param[in] name_file полное имя файла с расширением
 * @param[out] data массив std::vector
 * @return size_t - ::e_type_completion
 */
template <typename geo_elem>
size_t ReadGridGeo(const std::string& name_file, std::vector<geo_elem>& data) {
  FILE* f;
  OPEN_FILE(f, name_file.c_str(), "rb");

  int n;
  fread(&n, sizeof(int), 1, f);
  data.resize(n);
  for (auto& el : data) {
    fread(&el.geo, sizeof(el.geo), 1, f);
  }
  fclose(f);

  return e_completion_success;
}

/**
 * @brief Чтение нормалей к ячейкам сетки
 *
 * @param[in] name_file_normals полное имя файла со структурой нормалей
 * @note исходный файл генерируется функцией записи ::WriteNormals
 * @param[out] normals массив структур ::Normals
 * @return int ::e_type_completion
 */
int ReadNormals(const std::string& name_file_normals,
                std::vector<Normals>& normals);

/**
 * @brief Чтение данных, распределённых по ячейкам сетки
 *
 * @param[in] class_file_vtk конфигурация читаемых данных ::e_grid_vtk_config_t
 * @param[in] main_dir путь до папки с файлами данных
 * @param[out] density массив ::Type. плотность
 * @param[out] absorp_coef массив ::Type. коэффициент поглащения
 * @param[out] rad_en_loose_rate массив ::Type. коэффициент рассеяния
 * @param[out] velocity массив ::Vector3. скорость
 * @param[out] pressure массив ::Type. давление
 * @param[in] is_print флаг печати информации о считанных файлах
 * (default=false).
 * @note массивы могут не заполняться. В зависимости от class_file_vtk.
 * @return int ::e_type_completion
 */
int ReadData(const size_t class_file_vtk, const std::string& main_dir,
             std::vector<Type>& density, std::vector<Type>& absorp_coef,
             std::vector<Type>& rad_en_loose_rate,
             std::vector<Vector3>& velocity, std::vector<Type>& pressure,
             const bool is_print = false);

#ifdef ILLUM
/**
 * @brief Чтение данных, распределённых по ячейкам сетки
 * @warning НЕ ИСПОЛЬЗОВАТЬ!!!!
 * @param[in] mode конфигурация решателя
 * @param[in] main_dir путь до папки с файлами данных
 * @param[out] grid структура сетки
 * @return int ::e_type_completion
 */
int ReadData(const solve_mode_t& mode, const std::string& main_dir,
             grid_t& grid);

/**
 * @brief Чтение результатов трассировки.
 *
 * Чтение результатов трассировки модулей BUILD_GRAPH и BUILD_DATA_TO_ILLUM.
 * Учитывает конфигурацию MPI.
 *
 * @param[in] count_dir число направлений излучения
 * @param[in] gbl_files структура с конфигурацией входных файлов
 * @param[out] vec_x точки пересечения лучей с гранями тетраэдров
 * @param[out] face_states битовые флаги входящий/входящий граней
 * @param[out] vec_x0 точка начало луча в локальных координатах грани
 * @param[out] sorted_id_cell перенумерация ячеек вдоль луча
 * @param[out] vec_res_bound массив с предрасчитамм граничными значениями
 * @return int ::e_type_completion
 */
int ReadRadiationTrace(const int count_dir, const global_files_t& gbl_files,
                       std::vector<BasePointTetra>& vec_x,
                       std::vector<std::vector<int>>& face_states,
                       std::vector<std::vector<cell_local>>& vec_x0,
                       std::vector<std::vector<int>>& sorted_id_cell,
                       std::vector<Type>& vec_res_bound);
#endif  //! ILLUM

}  // namespace bin
}  // namespace files_sys

#endif  //! READER_BIN
