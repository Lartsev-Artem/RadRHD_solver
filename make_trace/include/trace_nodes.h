#if !defined TRACE_LOGIC_FUNC_H && defined MAKE_TRACE
#define TRACE_LOGIC_FUNC_H

#include "geo_types.h"
#include "trace_struct.h"

/**
 \todo в описание коротких характеристик
* \note алгоритм устроен так: для каждой выходящей грани надо знать три точки, чтобы построить на грани линейную функцию ax+by+c
 * функция находит координаты этих точек в глобальной системе координат/
*/

namespace trace {
int GetNodes(const int num_cell, const std::vector<Face> &grid, const ShortId num_cur_out_face,
             const Eigen::Matrix4d &vertex_tetra, const int face_state, const Vector3 &direction, const std::vector<Normals> &normals,
             const std::vector<int> &all_pairs_face, const BasePointTetra &x,
             std::vector<Type> &vec_res_bound, std::vector<cell_local> &vec_x0);

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