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

#include "aligned_vector.h"

extern solve_mode_t _solve_mode;
extern hllc_value_t _hllc_cfg;

/**
 * @brief Структура хранящая и организующая доступ к табулированной функция двух переменных
 *
 */
struct TableFunc {
  int Nx; ///< число узлов по первой переменной x
  int Ny; ///< число узлов по второй переменной y

  Type step_x; ///< шаг табуляции по x
  Type step_y; ///< шаг табуляции по y
  Type min_x;  ///< минимальное значение аргумента x
  Type min_y;  ///< минимальное значение аргумента y
  Type max_x;  ///< максимальное значение аргумента ч
  Type max_y;  ///< максимальное значение аргумента y

  std::vector<Type> data; ///< значения табулированной функции в формате (x*Ny + y)

  TableFunc(int nx = 0, int ny = 0) : Nx(nx), Ny(ny) {
    data.resize(nx * ny, 0);
  }

  /// \todo линейная интерполяция за пределы табуляции
  Type operator()(Type x, Type y); ///< безопасный доступ к значению функции
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

/**
 * @brief Геометрия грани
 *
 */
struct geo_face_t {

  uint32_t id_l : 25; ///< номер "левой" ячейки. всегда определена
  uint32_t id_l_node : 6;
  uint32_t is_regular : 1;

  int id_r : 26; ///< номер "правой" ячейки. включает признак границы
  uint32_t id_r_node : 6;

  Vector3 n; ///< нормаль к грани (не ориентирована)
  Type S;    ///< площадь грани

  geo_face_t() : id_l(0), id_r(0), n(Vector3(0, 0, 0)), S(0),
                 id_l_node(0), id_r_node(0), is_regular(1) {}
};

/**
 * @brief структура грани
 *
 */
struct face_t {
  flux_t f;       ///< поток определенный на грани
  geo_face_t geo; ///< геометрия ячейки
  face_t(){};
};

/**
 * @brief Геометрия ячейки
 *
 */
struct geo_cell_t {

  struct bits_flag {
    uint8_t bits;
    bits_flag(const uint8_t a = 0) : bits(a) {}
    bool operator[](int i) {
      return CHECK_BIT(bits, i);
    }
  };

  int id_faces[CELL_SIZE]; ///< номера связанных граней в нумерации face_t
  Type V;                  ///< объем ячейки
  Vector3 center;          ///< центр ячейки
  bits_flag sign_n;        ///< знак нормали соответствующей грани
  uint8_t node;

  geo_cell_t() : id_faces{0, 0, 0, 0}, V(0),
                 sign_n(0x7), node(0),
                 center(Vector3(0, 0, 0)) {}
};

/**
 * @brief структура определяющей ячейки для трассировки МКХ
 *
 */
struct cell_local // для каждой ячейки и каждого направления
{
#ifdef INTERPOLATION_ON_FACES
  Vector2 x0; ///< локальная координата входного узла для интерполяции
#endif
  Type s;             ///< расстояние x0-x
  ShortId in_face_id; ///< id выходной грани
};

/// \brief сжатая битовая структура для хранения локальных номеров граней
union face_loc_id_t {
  struct
  {
    uint8_t a : 2;
    uint8_t b : 2;
    uint8_t c : 2;
    uint8_t d : 2;
  };
  uint8_t bits;
};

struct align_cell_local // для каждой ячейки и каждого направления
{
#ifdef INTERPOLATION_ON_FACES
  Vector2 x0; ///< локальная координата входного узла для интерполяции
#endif
  std_ext::aligned_vector<Type, 32> s;   ///< расстояние x0-x
  std::vector<face_loc_id_t> in_face_id; ///< id выходной грани
};

/**
 * @brief  данные на сетке связанные с излучением
 *
 */
struct illum_value_t {
  Type absorp_coef;       ///< коэффициент поглощения (не ослабления!)
  Type scat_coef;         ///< коэффициент рассеяния
  Type rad_en_loose_rate; ///< источник равновесного излучения

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

  illum_value_t(const int num_dir = 0) : absorp_coef(0),
                                         rad_en_loose_rate(0)

#ifndef USE_CUDA
                                         ,
                                         energy(0),
                                         stream(Vector3::Zero()),
                                         impuls(Matrix3::Zero()),
                                         div_stream(0),
                                         div_impuls(Vector3::Zero()) {
    illum.resize(num_dir * CELL_SIZE, 0);
  }
#else
  {
  }
#endif
};

/**
 * @brief связанная пара ячейка-грань в графе обхода ячеек
 *
 */
union graph_pair_t {
  struct
  {
    uint32_t loc_face : 2; ///< локальный номер грани на которой проходит расчёт излучения
    IntId cell : 30;       ///< номер ячейки, сквозь которую идет трассировка
  };
  uint32_t bits;
};

struct full_phys_data_t {
  const flux_t *val;
  Type T;
  Type logT;
  Type lorenz;
  Type vel;
  Type alpha;
  Type betta;
#ifdef SPECTRUM
  Type cosf;
  void InitDirection(const Vector3 &dir);
#endif

  full_phys_data_t(const flux_t *src) { Init(src); }
  full_phys_data_t() {}
  void Init(const flux_t *src);
};

struct elem_t {
  flux_t phys_val;
  flux_t conv_val;

#if defined ILLUM
  full_phys_data_t *cell_data = nullptr;
  illum_value_t illum_val;
#endif

  geo_cell_t geo; //геометрия элемента

  // private:
  //	elem_t(const elem_t& el) {};
};

#ifndef USE_CUDA
struct grid_t {
  IdType size;
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
  void InitMemory(const IdType num_cells, const IdType num_directions);
};
#else
struct grid_t {
  IdType size;      ///< размер сетки (число ячеек)
  IdType size_face; ///< число граней
  IdType size_dir;  ///< число направлений

  std::vector<elem_t> cells;
  std::vector<face_t> faces;
#ifdef USE_MPI
  mpi_hllc_t *mpi_cfg;                              ///< конфиг mpi структуры
#endif

#ifndef TRANSFER_CELL_TO_FACE
  std::vector<std::vector<Vector3>> inter_coef_all; ///< коэффициенты интерполяции локальные для каждого потока
#else
#ifdef SPECTRUM
  std::vector<std::vector<std::vector<Type>>> inter_coef_all; ///< коэффициенты интерполяции локальные для каждого потока[thread[frq{cell}]]
#else
  std::vector<std::vector<Type>> inter_coef_all; ///< коэффициенты интерполяции локальные для каждого потока
#endif
#endif

#ifdef SPECTRUM
  IdType size_frq;
  std::vector<Type> spectrum;
#endif

#ifdef SEPARATE_GPU
  const int loc_illum_size = 10;
  int loc_illum_wptr = 0;
  int loc_illum_rptr = 0;
#endif

  Type *Illum;
  Type *scattering;

  IdType loc_size;
  IdType loc_shift;

  Type *divstream;
  Vector3 *divimpuls;

#ifdef ON_FULL_ILLUM_ARRAYS
  Type *energy;
  Vector3 *stream;
  Matrix3 *impuls;
#endif

#ifdef ILLUM
  //! \brief Инициализация памяти под структуру физических данных
  void InitFullPhysData();
#endif

  void InitMemory(const IdType num_cells, const grid_directions_t &dir_grid);

  grid_t() : size(0), size_face(0), size_dir(0), loc_size(0), loc_shift(0), Illum(nullptr), scattering(nullptr),
             divstream(nullptr), divimpuls(nullptr)
#ifdef ON_FULL_ILLUM_ARRAYS
             ,
             energy(nullptr), stream(nullptr), impuls(nullptr)
#endif
#ifdef SPECTRUM
             ,
             size_frq(0)
#endif
#ifdef USE_MPI
             ,
             mpi_cfg(nullptr)
#endif

  {
  }
  ~grid_t();
};

#endif // USE_CUDA
extern std::vector<bound_size_t> hllc_loc_size;

#endif //! SOLVERS_STRUCT_H