#if !defined SOLVERS_STRUCT_H
#define SOLVERS_STRUCT_H

#include "global_def.h"
#include "global_types.h"

extern solve_mode_t _solve_mode;

/**
 * @brief Структура потоков
 * @warning операторы написаны строго под hllc в той реализации, в которой есть
 * @todo переписать нормально. Решить вопрос с копированием
 */
struct flux_t {
  Type d;    ///< плотность
  Vector3 v; ///< скорость
  Type p;    ///< давление

  flux_t() : d(0), v(Vector3(0, 0, 0)), p(0) {}
  flux_t(const Type a, const Type b, const Type c, const Type d, const Type e) : d(a), v(Vector3(b, c, d)), p(e) {}

  flux_t operator+(const flux_t &x);
  flux_t operator+=(const flux_t &x);
  flux_t operator-=(const flux_t &x);
  flux_t operator*(const Type x);
  flux_t operator-(const flux_t &x);
  flux_t operator/(const Type x);

  // это временно для свзяи со старым кодом
  Type operator[](const int i) const;
  Type &operator[](const int i);
  Type operator()(const int i);

  // private:
  flux_t(const flux_t &f) : d(f.d), v(f.v), p(f.p) {}
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

struct illum_value_t {
  Type absorp_coef;
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
  grid_t() { size = 0; }
#endif
  void InitMemory(const uint32_t num_cells, const uint32_t num_directions);
};

extern std::vector<bound_size_t> hllc_loc_size;
#endif // USE_CUDA
#endif //! SOLVERS_STRUCT_H