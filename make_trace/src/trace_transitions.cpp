#ifdef MAKE_TRACE
#include "trace_transitions.h"

void trace::FromGlobalToLocalTetra(const Eigen::Matrix4d &vertex_tetra, const Vector3 &global_coord, Vector3 &local_coord) {
  // vertex_tetra -> [X;Y;Z;1]
  // возможно надо будет использовать transpose из-за инициализации матрицы перехода

  Eigen::Matrix4d vertex_tetra_inverse = vertex_tetra.inverse();
  // сразу с транспонированием
  for (int i = 0; i < 3; ++i) {
    local_coord[i] = 0;
    for (int j = 0; j < 3; ++j)
      local_coord[i] += vertex_tetra_inverse(i, j) * global_coord[j];
    local_coord[i] += vertex_tetra_inverse(i, 3);
  }

  // local_coord = vertex_tetra * global_coord;
}
void trace::FromLocalToGlobalTetra(const Eigen::Matrix4d &vertex_tetra, const Vector3 &local_coord, Vector3 &global_coord) {
  // vertex_tetra -> [X,Y,Z,1]

  Type eta4 = 1 - local_coord[0] - local_coord[1] - local_coord[2];
  for (int i = 0; i < 3; ++i) {
    global_coord[i] = 0;
    for (int j = 0; j < 3; ++j)
      global_coord[i] += vertex_tetra(i, j) * local_coord[j];
    global_coord[i] += vertex_tetra(i, 3) * eta4;
  }

  // написать в ручную т.к. преобразования известны, и 4я строка постоянна и меняться не должна
  // global_coord = vertex_tetra.inverse() * local_coord;
}

void trace::FromTetraToPlane(const Matrix3 &transform_matrix, const Vector3 &start_point, const Vector3 &tetra_coord, Vector3 &plane_coord) {
  plane_coord = transform_matrix * (tetra_coord - start_point);
}
void trace::FromPlaneToTetra(const Matrix3 &inverse_transform_matrix, const Vector3 &start_point, const Vector3 &plane_coord, Vector3 &tetra_coord) {
  tetra_coord = inverse_transform_matrix * plane_coord + start_point;
}

Vector3 GetInterpolationCoef(const Matrix3 &interpolation_nodes, const Vector3 &function_value) {
  return interpolation_nodes.partialPivLu().solve(function_value);
}
Vector3 GetInterpolationCoefInverse(const Matrix3 &interpolation_nodes, const Vector3 &function_value) {
  return interpolation_nodes * function_value;
}

#endif //! MAKE_TRACE