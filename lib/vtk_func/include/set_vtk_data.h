#if !defined SET_DATA_H && defined USE_VTK
#define SET_DATA_H

#include "dbgdef.h"
#include <typeinfo>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

/**
 * @brief Функция добавляет в структуру сетки поле данных
 *
 * @tparam Type тип данных
 * @param[in] name_data задаваемое название данных в сетки vtk
 * @param[in] data массив данных
 * @param[inout] u_grid сетка vtk
 * @return int  ::e_type_completion
 * @note записывает как скалярные так и векторные данные
 */
template <typename Type>
int SetDoubleVtkData(const std::string &name_data, const std::vector<Type> &data, vtkSmartPointer<vtkUnstructuredGrid> &u_grid) {

  vtkSmartPointer<vtkDoubleArray> data_vtk = vtkSmartPointer<vtkDoubleArray>::New();
  data_vtk->SetNumberOfComponents(sizeof(Type) / sizeof(double));
  data_vtk->SetName(name_data.c_str());

  switch (sizeof(Type) / sizeof(double)) {
  case 1:
    for (auto el : data) {
      data_vtk->InsertNextTuple((double *)&el);
    }
    u_grid->GetCellData()->AddArray(data_vtk);
    break;

  case 3:
    for (auto &el : data) {
      data_vtk->InsertNextTuple((double *)&el);
    }
    u_grid->GetCellData()->SetActiveVectors(name_data.c_str());
    u_grid->GetCellData()->SetVectors(data_vtk);

    break;

  case 9:
    for (auto &el : data) {
      data_vtk->InsertNextTuple((double *)&el);
    }
    u_grid->GetCellData()->SetActiveTensors(name_data.c_str());
    u_grid->GetCellData()->SetTensors(data_vtk);
    break;

  default:
    RETURN_ERR("unknow type data\n");
  }

  return e_completion_success;
}

/**
 * @brief Функция задаёт на сетки поля данных, из файлов решения
 *
 * @param[in] address_solution адрес с файлами решения
 * @param[inout] u_grid сетка
 * @return int ::e_type_completion
 */
int SetSolutionFromFileToVtk(const std::string &address_solution, vtkSmartPointer<vtkUnstructuredGrid> &u_grid);
#endif //! SET_DATA_H