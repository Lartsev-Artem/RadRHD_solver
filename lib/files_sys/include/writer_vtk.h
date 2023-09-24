/**
 * @file writer_vtk.h
 * @brief Файл содержит функции вывода vtk сеток
 * @warning требуется подключенная поддержка VTK
 *
 */

#ifndef WRITER_VTK
#define WRITER_VTK

#include "prj_config.h"

#ifdef USE_VTK

#include <vtkCellData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkGenericDataObjectWriter.h>
#include <vtkUnstructuredGrid.h>

/*! \addtogroup file_sys Файловый модуль
    @{
*/

namespace files_sys {

/**
 * @brief Пространство имён подмодуля vtk файлов
 *
 */
namespace vtk {

/**
 * @brief
 *
 * @tparam str_t строка
 * @param name_file_vtk полное имя файла с расширением vtk
 * @param grid сетка
 * @param format_bin флаг вывода в двоичном формате(default = true)
 * @return int ::e_type_completion
 */
template <typename str_t>
int WriteVtkGrid(str_t name_file_vtk, const vtkSmartPointer<vtkUnstructuredGrid> &grid, const bool format_bin = true) {

  vtkSmartPointer<vtkGenericDataObjectWriter> writer = vtkSmartPointer<vtkGenericDataObjectWriter>::New();

  if (format_bin) {
    writer->SetFileTypeToBinary();
  }

  writer->SetFileName(name_file_vtk.c_str());
  writer->SetInputData(grid);
  writer->Write();
  return 0;
}

} // namespace vtk
} // namespace files_sys

#endif //! USE_VTK
#endif // !WRITER_VTK