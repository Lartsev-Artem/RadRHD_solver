/**
 * @file trace_nodes.h
 * @brief Файл содержит функции поиска узлов интерполяции
 * @note Важно, что при переходе от глобального тетраэдра к локальному, вершины меняют свою нумерацию по правилу:
 * 0->3,
 * 1->2,
 * 2->0,
 * 3->1.
 *
 */
#if !defined TRACE_LOGIC_FUNC_H && defined MAKE_TRACE
#define TRACE_LOGIC_FUNC_H

#include "geo_types.h"
#include "solvers_struct.h"
#include "trace_struct.h"

namespace trace {

/**
 * @brief Функция находит три определяющих узла в локальных координатах для
 * трех узлов на выходящей гране текущей ячейки
 *
 * @param[in] num_cell номер ячейки
 * @param[in] num_out_face номер выходящей грани
 * @param[in] grid грани всех ячеек
 * @param[in] vertex_tetra матрица перехода в координаты локального тетраэдра
 * @param[in] face_state признаки входящих/выходящих граней
 * @param[in] direction направление
 * @param[in] normals структура нормалей к ячейке
 * @param[in] neighbours список соседей
 * @param[in] x узлы интерполяции тетраэдра в глобальных координатах
 * @param[out] vec_x0 список точек определяющих узлов в лок. системе координат
 * @return int ::e_type_completion
 * @note массив хранит упорядоченные данные (добавление через push_back())
 */
int GetLocNodes(const int num_cell, const ShortId num_out_face, const std::vector<Face> &grid,
                const Eigen::Matrix4d &vertex_tetra, const bits_flag_t face_state,
                const Vector3 &direction, const Normals &normals,
                const std::vector<IntId> &neighbours, const BasePointTetra &x,
                std::vector<cell_local> &vec_x0);
/**
 * @brief Функция возвращает точки интерполяции на всех гранях в глобальной системе координат
 *
 * @param[in] vertexs матрицы перехода в локальную систему тетраэдра
 * @param[out] vec_x координаты узлов интерполяции в глобальной системе координат
 * @return int ::e_type_completion
 */
int GetInterpolationNodes(const std::vector<Eigen::Matrix4d> &vertexs, std::vector<BasePointTetra> &vec_x);
} // namespace trace
#endif //! TRACE_LOGIC_FUNC_H