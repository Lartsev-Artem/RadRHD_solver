/**
 * @file writer_bin.h
 * @brief Запись бинарных данных

 */

#ifndef WRITER_BIN
#define WRITER_BIN

#include <string>
#include <vector>

#include "geo_types.h"
#include "solvers_struct.h"

/*! \addtogroup file_sys Файловый модуль
    @{
*/

namespace files_sys {
namespace bin {

// поэлементная запись данных в файл
#define WRITE_FILE_ELEM(name_file, data, value)  \
  {                                              \
    FILE *f;                                     \
    OPEN_FILE(f, name_file, "wb");               \
    int n = data.size();                         \
    fwrite(&n, sizeof(int), 1, f);               \
    for (auto &el : data) {                      \
      fwrite(&el.value, sizeof(el.value), 1, f); \
    }                                            \
    fclose(f);                                   \
  }

/**
 * @brief Запись бинарного файла
 *
 * @tparam Str_Type символьны тип
 * @tparam T тип записываемых данных
 * @param[in] name_file полное имя файла с расширением
 * @param[in] data массив std::vector
 * \note  Файл содержит в первой строке число элементов. Далее последовательные данные
 * @return int ::e_type_completion
 */
template <typename Str_Type, typename T>
int WriteSimple(const Str_Type name_file, const std::vector<T> &data) {

  FILE *f;
  OPEN_FILE(f, std::string(name_file).c_str(), "wb");

  int n = (int)data.size();
  fwrite(&n, sizeof(int), 1, f);
  fwrite(data.data(), sizeof(T), n, f);

  fclose(f);
  return e_completion_success;
}

/**
 * @brief Запись линейного массива размера  в файл
 *
 * @tparam Str_Type символьны тип
 * @tparam T тип записываемых данных
 * @param[in] name_file полное имя файла с расширением
 * @param[in] n размерность массива
 * @param[in] data массив T*
 * @return int  ::e_type_completion
 */
template <typename Str_Type, typename T>
int WriteSimple(const Str_Type name_file, const int n, const T *data) {

  if (data != nullptr) {
    FILE *f;
    OPEN_FILE(f, std::string(name_file).c_str(), "wb");

    fwrite(&n, sizeof(int), 1, f);
    fwrite(data, sizeof(T), n, f);
    fclose(f);
    return e_completion_success;
  }

  WRITE_LOG("no data for %s\n", std::string(name_file).c_str());
  return e_completion_fail;
}

/**
 * @brief Запись геометрии элементов сетки
 *
 * @tparam geo_elem - структура, содержащая поле "geo"
 * @param[in] name_file полное имя файла с расширением
 * @param[in] data массив std::vector
 * @return int - ::e_type_completion
 */
template <typename geo_elem>
int WriteGridGeo(const std::string &name_file, const std::vector<geo_elem> &data) {

  FILE *f;
  OPEN_FILE(f, name_file.c_str(), "wb");

  int n = (int)data.size();
  fwrite(&n, sizeof(int), 1, f);
  for (auto &el : data) {
    fwrite(&el.geo, sizeof(el.geo), 1, f);
  }
  fclose(f);

  return e_completion_success;
}

/**
 * @brief Запись структуры нормалей в файл
 *
 * @param[in] name_file_normals полное имя файла с расширением
 * @param[in] normals
 * @return int ::e_type_completion
 */
int WriteNormals(const std::string &name_file_normals, std::vector<Normals> &normals);

/**
 * @brief Запись файлов решения
 *
 * @param[in] main_dir путь к выходной директории
 * @param[in] grid сетка
 * @return int  ::e_type_completion
 */
int WriteSolution(const std::string &main_dir, const grid_t &grid);

#ifdef USE_MPI
#include "mpi_ext.h"
/**
 * @brief Запись файлов решения разделённого по узлам
 *
 * @param[in] main_dir путь к выходной директории
 * @param[in] grid сетка
 * @return int  ::e_type_completion
 */
int WriteFileSolutionMPI(const std::string &main_dir, const grid_t &grid);

/**
 * @brief Запись распределённого по узлам линейного массива размера n в файл
 *
 * @tparam Str_Type символьны тип
 * @tparam T тип записываемых данных
 * @param[in] file полное имя файла с расширением
 * @param[in] n размерность массива
 * @param[in] data массив T*
 * @param[in] left_offset левая граница данных на узле
 * @param[in] right_offset правая граница данных на узле
 * @return int  ::e_type_completion
 */
template <typename T, typename Tidx>
int WriteSimpleMPI(const std::string &file, const Tidx n, const T *data, const Tidx left_offset, const Tidx right_offset) {

  if (data != nullptr) {
    int myid = get_mpi_id();
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_File fh = MPI_FILE_NULL;
    MPI_File_open(comm, file.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

    DIE_IF(fh == MPI_FILE_NULL)

    if (myid == 0) {
      MPI_File_seek(fh, 0, MPI_SEEK_SET);
      MPI_File_write(fh, &n, 1, MPI_INT, MPI_STATUSES_IGNORE);
    } else {
      MPI_File_seek(fh, sizeof(int) + left_offset * sizeof(data[0]), MPI_SEEK_SET);
    }
    MPI_File_write(fh, data, (right_offset - left_offset) * sizeof(data[0]) / sizeof(double), MPI_DOUBLE, MPI_STATUSES_IGNORE);
    MPI_File_close(&fh);
    return e_completion_success;
  }
  WRITE_LOG("no data for %s\n", std::string(file).c_str());
  return e_completion_fail;
}

/**
 * @brief Запись данных распределённых по узлам
 *
 * @param[in] comm коммуникатор
 * @param[in] file_name имя файла
 * @param[in] cells сетка
 * @param[in] data_offset смещение данных в структуре elem_t
 * @param[in] size_elem размер записываемых данных
 * @param[in] left_offset левая граница данных на узле
 * @param[in] right_offset правая граница данных на узле
 * @warning Требует проверки
 */
template <typename type1>
void WriteFileVectorMPI(MPI_Comm comm, const std::string &file_name,
                        const std::vector<elem_t> &cells, uint32_t data_offset, uint32_t size_elem,
                        type1 left_offset, type1 right_offset) {
  int myid;
  MPI_Comm_rank(comm, &myid);

  MPI_File fh = MPI_FILE_NULL;
  MPI_File_open(comm, file_name.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
  if (fh == MPI_FILE_NULL)
    D_LD;

  if (myid == 0) {
    MPI_File_seek(fh, 0, MPI_SEEK_SET);
    int n = cells.size();
    MPI_File_write(fh, &n, 1, MPI_INT, MPI_STATUSES_IGNORE);
  } else {
    MPI_File_seek(fh, sizeof(int) + left_offset * size_elem, MPI_SEEK_SET);
  }

  for (int i = left_offset; i < right_offset; i++) {
    double *val = (double *)((uint8_t *)&cells[i] + data_offset);
    MPI_File_write(fh, val, size_elem / sizeof(double), MPI_DOUBLE, MPI_STATUSES_IGNORE);
  }
  MPI_File_close(&fh);
}

// поэлементная запись данных в файл с разных узлов
#define WRITE_FILE_ELEM_MPI(comm, file, cells, data, left_offset, right_offset)                                    \
  {                                                                                                                \
    MPI_File fh = MPI_FILE_NULL;                                                                                   \
    MPI_File_open(comm, file, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);                              \
    DIE_IF(fh == MPI_FILE_NULL)                                                                                    \
                                                                                                                   \
    if (myid == 0) {                                                                                               \
      MPI_File_seek(fh, 0, MPI_SEEK_SET);                                                                          \
      int n = cells.size();                                                                                        \
      MPI_File_write(fh, &n, 1, MPI_INT, MPI_STATUSES_IGNORE);                                                     \
    } else {                                                                                                       \
      MPI_File_seek(fh, sizeof(int) + left_offset * sizeof(cells[0].data), MPI_SEEK_SET);                          \
    }                                                                                                              \
                                                                                                                   \
    for (int i = left_offset; i < right_offset; i++) {                                                             \
      MPI_File_write(fh, &cells[i].data, sizeof(cells[0].data) / sizeof(double), MPI_DOUBLE, MPI_STATUSES_IGNORE); \
    }                                                                                                              \
    MPI_File_close(&fh);                                                                                           \
  }

#endif

} // namespace bin
} // namespace files_sys
#endif // !WRITER_BIN
