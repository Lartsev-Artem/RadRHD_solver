/**
 * @file reader_vtk.h
 * @brief Чтение vtk данных
 * @warning требуется подключенная поддержка VTK
 *
 */

#include "prj_config.h"

#if (!defined READER_VTK && defined USE_VTK)

#define READER_VTK

#include <vtkCellData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkGenericDataObjectWriter.h>
#include <vtkUnstructuredGrid.h>

#include "global_def.h"

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
 * @brief Чтение файла с сеткой в формате vtk
 *
 * @tparam vtk_grid тип сетки
 * @param[in] name_file полное имя файла с расширением
 * @param[out] grid сетка
 * @warning только тип сетки vtkUnstructuredGrid!
 * @return int ::e_type_completion
 */
template <typename vtk_grid>
int Read(const std::string &name_file, vtkSmartPointer<vtk_grid> &grid) {
  vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
      vtkSmartPointer<vtkGenericDataObjectReader>::New();
  reader_vtk->ReadAllScalarsOn();
  reader_vtk->SetFileName(name_file.c_str());
  reader_vtk->Update();

  if (reader_vtk->IsFileUnstructuredGrid()) {
    grid = reader_vtk->GetUnstructuredGridOutput();
    grid->Modified();
  } else {
    RETURN_ERR("Error read file_vtk\n file: %s is not UnstructuredGrid\n", name_file.c_str());
  }

  reader_vtk->GetOutput()->GlobalReleaseDataFlagOn();

  std::cout << "Grid has Cell: " << grid->GetNumberOfCells() << '\n';
  return e_completion_success;
}

/**
 * @brief
 *
 * @tparam vtk_grid тип сетки
 * @param[in] class_file_vtk конфигурация читаемых данных ::e_grid_vtk_config_t
 * @param[in] name_file_vtk полное имя файла с расширением vtk
 * @param[out] unstructured_grid сетка
 * @param[out] density массив плотности
 * @param[out] absorp_coef массив коэф. поглащения
 * @param[out] rad_en_loose_rate массив коэф. рассеяния
 * @param[in] is_print признак печати результатов чтения. default=false
 * @return int ::e_type_completion
 * @warning имена поля данных задаются в ручную!!!
 */
template <typename vtk_grid>
int Read(const size_t class_file_vtk, const std::string &name_file_vtk,
         vtkSmartPointer<vtk_grid> &unstructured_grid,
         vtkDataArray *&density, vtkDataArray *&absorp_coef,
         vtkDataArray *&rad_en_loose_rate, const bool is_print = false) {
  Read(name_file_vtk, unstructured_grid);

  switch (class_file_vtk) {
  case e_grid_cfg_default:
    density = NULL;
    absorp_coef = NULL;
    rad_en_loose_rate = NULL;
    break;

  case e_grid_cfg_radiation:
    density = unstructured_grid->GetCellData()->GetScalars("density");
    absorp_coef = unstructured_grid->GetCellData()->GetScalars("absorp_coef");
    rad_en_loose_rate =
        unstructured_grid->GetCellData()->GetScalars("radEnLooseRate");
    break;

  default:
    RETURN_ERR("Bad type vtk\n");
  }

  if (is_print) {
    std::cout << "Grid has " << unstructured_grid->GetNumberOfPoints()
              << " points.\n";
    std::cout << "Grid has " << unstructured_grid->GetNumberOfCells()
              << " cells.\n";
    if (class_file_vtk) {
      std::cout << "density_Size: " << density->GetSize() << '\n';
      std::cout << "absorp_coef_Size: " << absorp_coef->GetSize() << '\n';
      std::cout << "Q_Size: " << rad_en_loose_rate->GetSize() << '\n';
    }
  }

  return e_completion_success;
}
} // namespace vtk
} // namespace files_sys

#endif // USE_VTK && !READER_VTK
