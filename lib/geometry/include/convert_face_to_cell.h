
#ifndef CONVERT_FACE_TO_CELL_H
#define CONVERT_FACE_TO_CELL_H

#include "global_def.h"

/**
 * @brief Функция конвертирует значения с граней в значения на ячейке для выбранного направления
 *
 * @tparam Type тип данных с определёнными операторами /=, +=
 * @param[in] size_grid размер сетки
 * @param[in] num_dir номер направления
 * @param[in] illum_on_face массив скалярных данных на гранях
 * @param[in] zero ноль в смысле типа Type
 * @param[out] illum_in_cell массив на ячейке
 * @return int ::e_type_completion
 */
template <typename Type>
int GetDirectionDataFromFace(const int size_grid, const int num_dir, const Type *data_on_face, const Type zero, std::vector<Type> &data_in_cell) {

  if (data_on_face == nullptr)
    RETURN_ERR("data_on_face hasn't enough data\n");

  data_in_cell.assign(size_grid, zero);

  for (int j = 0; j < size_grid; j++) {
    const int N = num_dir * CELL_SIZE * size_grid + j * CELL_SIZE;
    for (int k = 0; k < CELL_SIZE; k++) {
      data_in_cell[j] += data_on_face[N + k];
    }
    data_in_cell[j] /= CELL_SIZE;
  }
  return e_completion_success;
}

/// \todo rename
template <typename Type>
int GetDirectionDataFromCellOrder(const int size_grid, const int size_dir, const int num_dir, const Type *data_on_face, const Type zero, std::vector<Type> &data_in_cell) {

  if (data_on_face == nullptr)
    RETURN_ERR("data_on_face hasn't enough data\n");

  data_in_cell.assign(size_grid, zero);

  for (int j = 0; j < size_grid; j++) {
    {
      data_in_cell[j] = data_on_face[j * size_dir + num_dir];
    }
  }
  return e_completion_success;
}

#endif //! CONVERT_FACE_TO_CELL_H