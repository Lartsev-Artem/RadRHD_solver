#ifndef GEO_TYPES
#define GEO_TYPES

#include <Eigen/Dense>

#include "global_def.h"

typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector2d Vector2;
typedef Eigen::VectorXd VectorX;
typedef Eigen::Matrix3d Matrix3;
typedef Eigen::Vector4d Vector4;
typedef Eigen::Matrix4d Matrix4;
typedef Eigen::MatrixXd MatrixX;
typedef double Type;
typedef uint8_t ShortId;
typedef int IntId;

///\todo: d в новый файл  с классами
struct Normals {
  std::vector<Vector3> n;
  Normals(const int size = CELL_SIZE) { n.resize(size); }
  ~Normals() { n.clear(); }
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
  int size;
  std::vector<direction_s> directions;
  Type full_area;
};

struct BasePointTetra //узлы интерполяции всех тетраэдров // В перспективе
                      //можно уйти к граням
{
  Vector3 x[CELL_SIZE][NUMBER_OF_MEASUREMENTS];

  Vector3 operator()(const int i, const int j) const { return x[i][j]; }
};
struct cell_local // для каждой ячейки и каждого направления
{
  Type s;             //расстояние x0-x
  Vector2 x0;         //локальная координата входного узла для интерполяции
  ShortId in_face_id; // id выходной грани

  friend std::ostream &operator<<(std::ostream &out, const cell_local &point);
};

#include "solvers_types.h"
struct grid_t {
  struct elem_t {
    flux_t phys_val;
  };

  struct face_t {
    int a;
  };

  int size;
  std::vector<elem_t> cells;
  std::vector<face_t> faces;

  /// \todo all config!
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

#endif // GEO_TYPES