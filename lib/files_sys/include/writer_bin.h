/**
 * @file writer_bin.h
 * @brief Запись бинарных данных

 */

#ifndef WRITER_BIN
#define WRITER_BIN

#include <string>
#include <vector>

#include "geo_types.h"
#include "solvers_types.h"

/*! \addtogroup file_sys Файловый модуль
    @{
*/

namespace files_sys {
namespace bin {

/**
 * @brief Запись бинарного файла
 *
 * @tparam Str_Type символьны тип
 * @tparam T тип записываемых данных
 * @param[in] name_file полное имя файла с расширением
 * @param[in] data массив std::vector
 * \note  Файл содержит в первой строке число элементов. Далее последовательные данные
 * @return size_t ::e_type_completion
 */
template <typename Str_Type, typename T>
size_t WriteSimple(const Str_Type name_file, const std::vector<T> &data) {

  FILE *f;
  OPEN_FILE(f, std::string(name_file).c_str(), "wb");

  int n = data.size();
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
 * @return size_t
 */
template <typename Str_Type, typename T>
size_t WriteSimple(const Str_Type name_file, const int n, const T *data) {

  FILE *f;
  OPEN_FILE(f, std::string(name_file).c_str(), "wb");

  fwrite(&n, sizeof(int), 1, f);
  fwrite(data, sizeof(T), n, f);
  fclose(f);
  return e_completion_success;
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

  int n = data.size();
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
 * @return size_t  ::e_type_completion
 */
size_t WriteSolution(const std::string &main_dir, const grid_t &grid);

/// \todo CHECK THIS!!!
#ifdef RHLLC_MPI
#error "todo module"
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

#define WRITE_FILE_VECTOR_MPI(comm, file, cells, data, left_offset, right_offset)                                  \
  {                                                                                                                \
    MPI_File fh = MPI_FILE_NULL;                                                                                   \
    MPI_File_open(comm, file, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);                              \
    if (fh == MPI_FILE_NULL)                                                                                       \
      D_LD;                                                                                                        \
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

#define WRITE_FILE_MPI(comm, file, data, n, left_offset, right_offset)                                                          \
  {                                                                                                                             \
    MPI_File fh = MPI_FILE_NULL;                                                                                                \
    MPI_File_open(comm, file, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);                                           \
    if (fh == MPI_FILE_NULL)                                                                                                    \
      D_LD;                                                                                                                     \
                                                                                                                                \
    if (myid == 0) {                                                                                                            \
      MPI_File_seek(fh, 0, MPI_SEEK_SET);                                                                                       \
      MPI_File_write(fh, &n, 1, MPI_INT, MPI_STATUSES_IGNORE);                                                                  \
    } else {                                                                                                                    \
      MPI_File_seek(fh, sizeof(int) + left_offset * sizeof(data[0]), MPI_SEEK_SET);                                             \
    }                                                                                                                           \
    MPI_File_write(fh, data, (right_offset - left_offset) * sizeof(data[0]) / sizeof(double), MPI_DOUBLE, MPI_STATUSES_IGNORE); \
    MPI_File_close(&fh);                                                                                                        \
  }
#endif

} // namespace bin
} // namespace files_sys
#endif // !WRITER_BIN
