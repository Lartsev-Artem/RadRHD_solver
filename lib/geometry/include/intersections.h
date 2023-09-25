#ifndef INTERSECTIONS_H
#define INTERSECTIONS_H
#include "geo_types.h"
#include "global_types.h"

/**
 * @brief Пространство имён геометрических пересечений
 *
 */
namespace intersection {

void IntersectionWithPlane(const Face &face, const Vector3 &start_point, const Vector3 &direction, Vector3 &result);
int InTriangle(int number_face, const Face &cell_face, const Normals &normals_cell, const Vector3 &XX);

/**
 * @brief Определяет типы граней ячейки по отношению к направлению
 *
 * @param direction[in] направление
 * @param normals_cell[in] нормали к ячейке
 * @param face_type[out] тип ячейки ::e_face_in_out_type_t
 */
void FindInAndOutFaces(const Vector3 &direction, const Normals &normals_cell, bits_flag_t &face_type);
} // namespace intersection
#endif //! INTERSECTIONS_H