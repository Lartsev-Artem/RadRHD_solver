/**
 * @file geo_types.h
 * @brief Объявления геометрических структур и типов
 *
 */
#ifndef GEO_TYPES
#define GEO_TYPES

// #if !__NVCC__
// #define EIGEN_NO_CUDA
// #endif
// workaround issue between gcc >= 4.7 and cuda 5.5 (совместимость с компиляторами (см. док. Eigen3: https://eigen.tuxfamily.org/dox/TopicCUDA.html))
#if (defined __GNUC__) && (__GNUC__ > 4 || __GNUC_MINOR__ >= 7)
#undef _GLIBCXX_ATOMIC_BUILTINS
#undef _GLIBCXX_USE_INT128
#endif
#include <Eigen/Dense>

#include "global_def.h"
#include "global_types.h"

#include <array>

typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector2d Vector2;
typedef Eigen::VectorXd VectorX;
typedef Eigen::Matrix3d Matrix3;
typedef Eigen::Vector4d Vector4;
typedef Eigen::Matrix4d Matrix4;
typedef Eigen::MatrixXd MatrixX;

struct Ray_t {
  Vector3 orig;
  Vector3 direction;

  Ray_t() {}
  Ray_t(const Vector3 o, const Vector3 dir) : direction(dir), orig(o) {}
};

struct Intersection_t {
  int id;
  Type dist;
  Vector3 point;
  Intersection_t(const int i = -1, const Type d = -1, const Vector3 &p = Vector3::Zero()) : id(i), dist(d), point(p) {}
};

///\todo: d в новый файл  с классами.
struct Normals {
  std::array<Vector3, CELL_SIZE> n;
};

struct Face {
  Vector3 A;
  Vector3 B;
  Vector3 C;
  Face &operator=(const Face &face) {
    A = face.A;
    B = face.B;
    C = face.C;
    return *this;
  }

  Vector3 &operator[](const int i) {
    return *((Vector3 *)((uint8_t *)&A + sizeof(Vector3) * i));
  }

  Face() {}
  Face(const Face &f) : A(f.A), B(f.B), C(f.C) {}
};

struct FaceCell {
  int face_id;
  Face face;
  FaceCell(const int id = 0, const Face &face_init = Face())
      : face_id(id), face(face_init) {}
};
std::ostream &operator<<(std::ostream &os, const std::pair<const int, FaceCell> &f); ///< ключ не выводим
std::ostream &operator<<(std::ostream &os, const FaceCell &f);
std::istream &operator>>(std::istream &is, FaceCell &f);

struct direction_s {
  Vector3 dir;
  Type area;
};

struct grid_directions_t {
  IdType size;
  IdType loc_size;
  IdType loc_shift;
  std::vector<direction_s> directions;
  Type full_area;
  grid_directions_t(const IdType N = 0) : size(N), loc_size(N), loc_shift(0), directions(N) {}
};

//это узлы интерполяции на гранях
struct BasePointTetra //узлы интерполяции всех тетраэдров // В перспективе можно уйти к граням
{
  Vector3 x[CELL_SIZE][NUMBER_OF_MEASUREMENTS];

  Vector3 operator()(const int i, const int j) const { return x[i][j]; }
};
/**
 * @brief Тип соседней ячейки к грани текущей
 * @note величины от 0 до N обозначают номер соседнего элемента
 *
 */
enum e_neighbor_code_t {
  e_neigh_code_undef = -5,     ///< сосед не определён
  e_neigh_code_in_bound = -2,  ///< внутрення граница области
  e_neigh_code_out_bound = -1, ///< внешняя граница
};

/**
 * @brief тип грани ячейки по отношению к направлению
 *
 */
enum e_face_in_out_type_t {
  e_face_type_out = 0, /// <выходящая грань
  e_face_type_in = 1   ///< входящая грань
};

/**
 * @brief Структура описывающая стандартный тетраэдр с
 * единичными размерами и прямыми углами и началом в точке {0,0,0}
 *
 */
struct BaseTetra_t {
  Vector3 start_point_plane_coord;  ///< начало координат плоскости
  Matrix3 transform_matrix;         ///< матрица перехода из базового тетраэдра в плоскость
  Matrix3 inverse_transform_matrix; ///< матрица перехода из плоскости в базовый тетраэдр

  Matrix3 straight_face; ///< 3 узла интерполяции
  Matrix3 inclined_face; ///< 3 узла интерполяции на наклонной плоскости

  Matrix3 straight_face_inverse; ///< 3 узла интерполяции
  Matrix3 inclined_face_inverse; ///< 3 узла интерполяции на наклонной плоскости

  //	Vector3 center_local_sphere;  // центр описанной сферы около стандартного тетраэдра

  BaseTetra_t() {
    // 3 узла интерполяции
    {
      straight_face << 1. / 6, 1. / 6, 1,
          2. / 3, 1. / 6, 1,
          1. / 6, 2. / 3, 1;
    }

    // 3 узла интерполяции на наклонной плоскости
    {
      inclined_face << 0, sqrt(2. / 3), 1,
          sqrt(2) / 4, 1. / (2 * sqrt(6)), 1,
          -sqrt(2) / 4, 1. / (2 * sqrt(6)), 1;
    }

    //Матрицы перехода из стандартного тетраэдра в координаты наклонной плоскости
    {
      transform_matrix << -1. / sqrt(2), 1. / sqrt(2), 0,
          -1. / sqrt(6), -1. / sqrt(6), sqrt(2. / 3),
          1. / sqrt(3), 1. / sqrt(3), 1. / sqrt(3);
    }

    //Матрицы перехода из наклонной плоскости в  координаты стандартного тетраэдра
    {
      inverse_transform_matrix << -1. / sqrt(2), -1. / sqrt(6), 1. / sqrt(3),
          1. / sqrt(2), -1. / sqrt(6), 1. / sqrt(3),
          0, sqrt(2. / 3), 1. / sqrt(3);
    }

    // Начало координата плоскости
    start_point_plane_coord << 0.5, 0.5, 0;

    inclined_face_inverse = inclined_face.inverse(); // в решении
    straight_face_inverse = straight_face.inverse(); // в решении
  }
};

#endif // GEO_TYPES