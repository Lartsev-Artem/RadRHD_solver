/**
 * @file solvers_struct.h
 * @brief Структуры связанные непосредственно с решателями излучения и газа
 *
 */
#if !defined SOLVERS_STRUCT_H
#define SOLVERS_STRUCT_H

#include "geo_types.h"
#include "global_def.h"
#include "global_types.h"

#include "solvers_config.h"

extern solve_mode_t _solve_mode;
extern hllc_value_t _hllc_cfg;

struct TableFunc {
  int Nx, Ny;

  Type step_x, step_y;
  Type min_x, min_y;
  Type max_x, max_y;

  std::vector<Type> data;

  TableFunc(int nx = 0, int ny = 0) : Nx(nx), Ny(ny) {
    data.resize(nx * ny, 0);
  }

  Type operator()(Type x, Type y);
};

/**
 * @brief Структура потоков
 * @warning операторы написаны строго под hllc в той реализации, в которой есть
 */
struct flux_t {
  Type d;    ///< плотность
  Vector3 v; ///< скорость
  Type p;    ///< давление

  flux_t() : d(0), v(Vector3(0, 0, 0)), p(0) {}
  flux_t(const Type a, const Vector3 b, const Type c) : d(a), v(b), p(c) {}
  flux_t(const flux_t &f) : d(f.d), v(f.v), p(f.p) {}

  flux_t operator+(const flux_t &x) const;
  void operator=(const flux_t &x);
  void operator+=(const flux_t &x);
  void operator-=(const flux_t &x);
  void operator*=(const Type x);
  void operator/=(const Type x);

  Type operator[](const int i) const;
  Type &operator[](const int i);
};

struct flux_all_t {
  flux_t phys_val;
  flux_t conv_val;
};
struct bound_size_t {
  int left;
  int right;
};

struct geo_face_t {
  int id_l;
  int id_r;

  Vector3 n;
  Type S;

  geo_face_t() : id_l(0), id_r(0), n(Vector3(0, 0, 0)), S(0) {}
};

struct face_t {
  flux_t f;
  geo_face_t geo; // геометрия ячейки
  face_t(){};
};

struct geo_cell_t {
  int id_faces[CELL_SIZE];
  Type V;
  bool sign_n[CELL_SIZE];
  Vector3 center;

  geo_cell_t() : id_faces{0, 0, 0, 0}, V(0),
                 sign_n{true, true, true, true},
                 center(Vector3(0, 0, 0)) {}
};

struct cell_local // для каждой ячейки и каждого направления
{
#ifdef INTERPOLATION_ON_FACES
  Vector2 x0; ///< локальная координата входного узла для интерполяции
#endif
  Type s;             ///< расстояние x0-x
  ShortId in_face_id; ///< id выходной грани
};

struct illum_value_t {
  Type absorp_coef;
  Type scat_coef;
  Type rad_en_loose_rate;

//в противном случае эти структуры вынесены как указатели в сетку и не привязаны к ячейкам
#if !defined USE_CUDA
  std::vector<Type> illum; // num_dir*CELL_SIZE

  Type energy;
  Vector3 stream;

  //! в перспективе не для квазистационарного уравнения
  // Type prev_energy;
  // Vector3 prev_stream;

  Matrix3 impuls;
  Type div_stream;
  Vector3 div_impuls;
#endif

  illum_value_t(const int num_dir = 0);
};

struct elem_t {
  flux_t phys_val;
  flux_t conv_val;

#if defined ILLUM
  illum_value_t illum_val;
#endif

  geo_cell_t geo; //геометрия элемента

  // private:
  //	elem_t(const elem_t& el) {};
};

#ifndef USE_CUDA
struct grid_t {
  int size;
  std::vector<elem_t> cells;
  std::vector<face_t> faces;

#if defined ILLUM
  Type *Illum;
  Type *scattering;
  std::vector<std::vector<Vector3>> inter_coef_all; ///< коэффициенты интерполяции локальные для каждого потока

  grid_t() : size(0), Illum(nullptr), scattering(nullptr) {}
  ~grid_t();

#else
  grid_t() : size(0) {}
#endif
  void InitMemory(const uint32_t num_cells, const uint32_t num_directions);
};
#else
struct grid_t {
  int size;
  std::vector<elem_t> cells;
  std::vector<face_t> faces;
  std::vector<std::vector<Vector3>> inter_coef_all; ///< коэффициенты интерполяции локальные для каждого потока

  Type *Illum;
  Type *scattering;

  int loc_size;
  int loc_shift;

  Type *divstream;
  Vector3 *divimpuls;

#ifdef ON_FULL_ILLUM_ARRAYS
  Type *energy;
  Vector3 *stream;
  Matrix3 *impuls;
#endif

  void InitMemory(const uint32_t num_cells, const uint32_t num_directions);

  grid_t() : size(0), loc_size(0), loc_shift(0), Illum(nullptr), scattering(nullptr),
             divstream(nullptr), divimpuls(nullptr)
#ifdef ON_FULL_ILLUM_ARRAYS
             ,
             energy(nullptr), stream(nullptr), impuls(nullptr)
#endif
  {
  }
  ~grid_t();
};

#endif // USE_CUDA
extern std::vector<bound_size_t> hllc_loc_size;

#endif //! SOLVERS_STRUCT_H