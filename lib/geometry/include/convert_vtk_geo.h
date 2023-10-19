/**
 * @file convert_vtk_geo.h
 * @brief Файл содержит методы пересчёта vtk формата во внутренний формат
 *
 */

#ifndef CONVERT_VTK_GEO_H
#define CONVERT_VTK_GEO_H

#include "prj_config.h"
#ifdef USE_VTK

#include <map>
#include <set>
#include <vector>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "geo_types.h"

/**
 * @brief Get the Faces Points object
 *
 * @param[in] unstructured_grid сетка в формате vtk
 * @param[out] faces массив с гранями по всем ячейкам (в глобальной нумерации (size_cell*number_of_cell + face))
 * @return int ::e_type_completion
 */
int GetFacesPoints(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Face> &faces);

/**
 * @brief Функция переводит заданные грани из vtk формата во внутренний
 *
 * @tparam stdT
 * @param[in] id_faces std контейнер id ячеек для перевода (в формате id_cell*cell_size + loc_id_face)
 * @param[in] unstructured_grid сетка в формате vtk
 * @param[out] faces контейнер с перестроенными гранями
 * @return int ::e_type_completion
 */
template <typename stdT>
int GetFacesById(const stdT &id_faces, const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::map<IntId, FaceCell> &faces) {

  faces.clear();

  FaceCell face_buf;
  vtkPoints *points_face;

  for (auto id : id_faces) {

    points_face = unstructured_grid->GetCell(id / CELL_SIZE)->GetFace(id % CELL_SIZE)->GetPoints();

    face_buf.face_id = id;
    for (int i = 0; i < CELL_SIZE - 1; i++) {
      face_buf.face[i] = Vector3(points_face->GetPoint(i));
    }

    faces.emplace(id, face_buf);
  }

  return e_completion_success;
}

/**
 * @brief Возвращает уникальный набор номеров граничных ячеек
 *
 * @param[in] unstructured_grid сетка в формате vtk
 * @param[out] boundary_cells набор(std::set) граничных ячеек
 * @return int ::e_type_completion
 */
int GetBoundaryCells(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::set<IntId> &boundary_cells);

/**
 * @brief Возвращает уникальный набор номеров граничных граней внутренней сферы
 *
 * @param[in] unstructured_grid сетка в формате vtk
 * @param[out] inter_boundary_faces набор(std::set) граничных граней
 * @return int ::e_type_completion
 */
int GetInterBoundaryFacesOfSphere(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::set<IntId> &inter_boundary_faces);

/**
 * @brief Возвращает  номеров граничных граней (с включением внутренней границы)
 * в глобальной нумерации (cell*cell_size + face)
 *
 * @param[in] unstructured_grid сетка в формате vtk
 * @param[out] boundary_faces номера граничных граней
 * @return int ::e_type_completion
 */
int GetBoundaryFacesId(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<IntId> &boundary_faces);
// #if NUMBER_OF_MEASUREMENTS == 3

/**
 * @brief Функция находит всех соседей 3d сетки (без привязки к геометрии)
 *
 * @param[in] unstructured_grid сетка в формате vtk
 * @param[out] neighbors массив соседей
 * @return int ::e_type_completion
 */
int GetNeighborFace3D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<int> &neighbors);

/**
 * @brief Функция определяет внешнии нормали и площади граней 3d сетки
 *
 * @param[in] unstructured_grid  сетка в формате vtk
 * @param[out] normals массив нормалей
 * @param[out] areas массив площадей
 * @return int ::e_type_completion
 */
int GetNormalAndAreas3D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Normals> &normals, std::vector<Type> &areas);

/**
 * @brief Функция определяет объёмы ячеек 3d сетки
 *
 * @param[in] unstructured_grid сетка в формате vtk
 * @param[out] volumes массив объёмов
 * @return int ::e_type_completion
 */
int GetVolume3D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Type> &volumes);

/**
 * @brief  Функция определяет центры ячеек 3d сетки
 *
 * \note под центром понимается центр вписанной сферы
 * @param[in] unstructured_grid сетка в формате vtk
 * @param[out] centers массив точек-центров
 * @return int ::e_type_completion
 */
int GetCentersOfCells3D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Vector3> &centers);

/**
 * @brief Функция определяет центры граней 3d сетки
 *
 * \note под центром понимается центр вписанной окружности
 * @param unstructured_grid сетка в формате vtk
 * @param centers массив точек-центров (в пространственных координатах)
 * @return int ::e_type_completion
 */
int GetCentersOfFaces3D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Vector3> &centers);

/**
 * @brief Функция формирует матрицу с координатами вершин ячейки и единицами в последней строке
 *
 * \note 4 вершины тетраэдра(по столбцам и единицы в нижний строке)
 *  применяется для перехода в локальные координаты грани при трассировке луча
 * @param[in] number_cell номер ячейки
 * @param[in] unstructured_grid сетка в формате vtk
 * @param[out] vertex_tetra матрица
 * @return int ::e_type_completion
 */
int GetVertexMatrix(const size_t number_cell, const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, Eigen::Matrix4d &vertex_tetra);
// #else

/**
 * @brief Функция находит всех соседей 2d сетки (без привязки к геометрии)
 *
 * @param[in] unstructured_grid сетка в формате vtk
 * @param[out] neighbors массив соседей
 * @return int::e_type_completion
 */
int GetNeighborFace2D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<int> &neighbors);

/**
 * @brief Функция определяет внешнии нормали и площади граней 2d сетки
 *
 * @param[in] unstructured_grid  сетка в формате vtk
 * @param[out] normals массив нормалей
 * @param[out] areas массив площадей
 * @return int ::e_type_completion
 */
int GetNormalAndAreas2D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Normals> &normals, std::vector<Type> &areas);

/**
 * @brief Функция определяет объёмы ячеек 2d сетки
 *
 * @param[in] unstructured_grid
 * @param[out] volumes
 * @return int ::e_type_completion
 */
int GetVolume2D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Type> &volumes);

/**
 * @brief Функция определяет центры ячеек 2d сетки
 *
 * \note под центром понимается центр вписанной окружности
 * @param[in] unstructured_grid сетка в формате vtk
 * @param[out] centers массив точек-центров
 * @return int ::e_type_completion
 */
int GetCentersOfCells2D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Vector3> &centers);

/**
 * @brief Функция возвращает ячейки сетки - поверхности
 *
 * @param[in] unstructured_grid сетка поверхности в формате vtk
 * @param[out] faces ячейки на поверхности
 * @return int ::e_type_completion
 */
int GetCellsPointsSurface(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Face> &faces);
// #endif // NUMBER_OF_MEASUREMENTS

#endif // USE_VTK
#endif //! CONVERT_VTK_GEO_H