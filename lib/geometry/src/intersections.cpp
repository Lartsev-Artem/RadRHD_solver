#include "intersections.h"

/**
 * @brief Построение локального базиса плоскости
 *
 * @details по начальной точке и нормали строит локальный базис картинной плоскости (vec1, vec2).
    нормаль дана.
    задаем второй вектор ортогонально, по св-ву скалярного умножения
    третий вектор из векторного произведения
 *
 * @param[in] normal нормаль к плоскости
 * @param[out] basis новый базис (тройка векторов)
 */
static void SetBasis(const Vector3 &normal, Matrix3 &basis) {
  Vector3 vec_1;
  Vector3 vec_2;

  if (fabs(normal[1]) < 1e-20) {
    vec_1[0] = 0;
    vec_1[1] = 1;
    vec_1[2] = 0;
  } else {
    vec_1[0] = 1;
    vec_1[2] = 0;
    vec_1[1] = -(normal[0] * vec_1[0] + normal[2] * vec_1[2]) / normal[1]; //св-во скалярного произведения (N, vec1)==0
  }

  // правильная ориентация базиса плоскости
  if (normal[1] < 0)
    vec_1 = -vec_1;

  // обычное векторное умножение.
  vec_2 = -(normal.cross(vec_1));

  vec_1.normalize();
  vec_2.normalize();

  basis.row(0) = vec_1;
  basis.row(1) = vec_2;
  basis.row(2) = normal;
}

/**
 * @brief Переход в локальный базис {start, local_basis}
 *
 * @param[in] start начало новой системы координат
 * @param[in] local_basis тройка векторов нового базиса
 * @param[in] point точка в 3d
 * @param[out] new_point точка в локальных координатах
 */
static void TransformTo2d(const Vector3 &start, const Matrix3 &local_basis,
                          const Vector3 &point, Vector3 &new_point) {

  new_point = Vector3::Zero();

  //перевод 3d точки в 2d (в локальном базисе {start, local_basis})
  for (int k = 0; k < 3; k++) {
    new_point[0] += (point[k] - start[k]) * local_basis(0, k);
    new_point[1] += (point[k] - start[k]) * local_basis(1, k);
  }
}

int intersection::IntersectionWithPlaneDisk(const Vector3 &X0, const Vector3 &n, Vector3 &res) {

  //  ----------полный расчет. Т.к. диск задается постоянной плоскостью, параметры можно задатб явно--------------
  /*
   {
  std::vector<Vector3> curface(3);		 // точки задающие плоскость диска
                  curface[0][0] = 1;
                  curface[0][1] = 0;
                  curface[0][2] = 0;

                  curface[1][0] = 0;//0;
                  curface[1][1] = 0.9928768384869221;//
                  curface[1][2] = 0.11914522061843064;//;

                  curface[2][0] = 2;//;
                  curface[2][1] = 0;
                  curface[2][2] = 0;// ;   // Wolfram
                  }

  * const std::vector<Vector3>& face,
  Vector3 A = face[0];
  Vector3 B = face[1];
  Vector3 C = face[2];

  Type a = A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2]);
  Type b = A[0] * (C[2] - B[2]) + B[0] * (A[2] - C[2]) + C[0] * (B[2] - A[2]);
  Type c = A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]);
  Type d = A[0] * (C[1] * B[2] - B[1] * C[2]) + B[0] * (A[1] * C[2] - C[1] * A[2]) + C[0] * (B[1] * A[2] - A[1] * B[2]);

  Type t = -(a * X0[0] + b * X0[1] + c * X0[2] + d) / (a * n[0] + b * n[1] + c * n[2]);
  */

  /*
  a= 0
  b= 0.1191452206184306
  c= -0.9928768384869221
  d= 0
  */

  const Type b = 0.1191452206184306;
  const Type c = -0.9928768384869221;

  const Type t = -(b * X0[1] + c * X0[2]) / (b * n[1] + c * n[2]);

  res = (t * n + X0);

  /*for (size_t i = 0; i < 3; i++)
          res[i] = (n[i] * t + X0[i]);*/

  return 0;
}

void intersection::IntersectionWithPlane(const Face &face, const Vector3 &start_point, const Vector3 &direction, Vector3 &result) {

  Type a, b, c, d; // параметры уравнения плоскости
  Type t;

  a = face.A[1] * (face.B[2] - face.C[2]) + face.B[1] * (face.C[2] - face.A[2]) + face.C[1] * (face.A[2] - face.B[2]);
  b = face.A[0] * (face.C[2] - face.B[2]) + face.B[0] * (face.A[2] - face.C[2]) + face.C[0] * (face.B[2] - face.A[2]);
  c = face.A[0] * (face.B[1] - face.C[1]) + face.B[0] * (face.C[1] - face.A[1]) + face.C[0] * (face.A[1] - face.B[1]);
  d = face.A[0] * (face.C[1] * face.B[2] - face.B[1] * face.C[2]) + face.B[0] * (face.A[1] * face.C[2] - face.C[1] * face.A[2]) + face.C[0] * (face.B[1] * face.A[2] - face.A[1] * face.B[2]);

  t = -(a * start_point[0] + b * start_point[1] + c * start_point[2] + d) / (a * direction[0] + b * direction[1] + c * direction[2]);

  for (int i = 0; i < 3; ++i)
    result[i] = (direction[i] * t + start_point[i]); // точка пересечения луча  (start->direction) с плоскостью!!! face
}

int intersection::InTriangle(int number_face, const Face &cell_face, const Normals &normals_cell, const Vector3 &XX) {
  /*face --- треугольник, X --- точка для проверки*/

  Vector3 A, B, C, X; // новые точки на плоскости
  {
    Eigen::Matrix3d basis;
    Vector3 n = normals_cell.n[number_face % 4];
    SetBasis(n, basis);
    TransformTo2d(cell_face.A, basis, cell_face.A, A);
    TransformTo2d(cell_face.A, basis, cell_face.B, B);
    TransformTo2d(cell_face.A, basis, cell_face.C, C);
    TransformTo2d(cell_face.A, basis, XX, X);
  }

  // линейная алгебра
  Type r1 = (A[0] - X[0]) * (B[1] - A[1]) - (B[0] - A[0]) * (A[1] - X[1]);
  Type r2 = (B[0] - X[0]) * (C[1] - B[1]) - (C[0] - B[0]) * (B[1] - X[1]);
  Type r3 = (C[0] - X[0]) * (A[1] - C[1]) - (A[0] - C[0]) * (C[1] - X[1]);

  if (r1 < 0 && r2 < 0 && r3 < 0)
    return true;
  else if (r1 > 0 && r2 > 0 && r3 > 0)
    return true;
  else
    return false;
}

void intersection::FindInAndOutFaces(const Vector3 &direction, const Normals &normals_cell, bits_flag_t &face_type) {
  face_type = 0;
  for (int i = 0; i < CELL_SIZE; ++i) {
    if (normals_cell.n[i].dot(direction) < -1e-10) {
      face_type = SET_BIT(face_type, i);
    }
  }
}
