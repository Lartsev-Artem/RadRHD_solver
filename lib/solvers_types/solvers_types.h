#ifndef SOLVERS_TYPES
#define SOLVERS_TYPES

#include "global_def.h"
#include "global_types.h"
#include "solvers_config.h"

struct flux_t {
public:
  Type d;
  Vector3 v;
  Type p;

  flux_t();
  flux_t(const Type a, const Type b, const Type c, const Type dd, const Type e);

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
  flux_t(const flux_t &f);
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

  geo_face_t();
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

  geo_cell_t();
};

struct illum_value_t {
  // Type density;
  Type absorp_coef;
  Type rad_en_loose_rate;

#if !defined USE_CUDA
  std::vector<Type> illum; // num_dir*CELL_SIZE

  // std::vector<Type> int_scattering; // (count_cells* count_directions, 0);
  Type energy;
  Vector3 stream;

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

extern std::vector<bound_size_t> hllc_loc_size;
#endif //! SOLVERS_TYPES