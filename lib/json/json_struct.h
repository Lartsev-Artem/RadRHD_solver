#ifndef JSON_STRUCT
#define JSON_STRUCT

#include "json.hpp"
typedef double hllc_value_type;

#define HLLC_VALUE_T                                                                                     \
  X(hllc_value_type, T) /*временной интервал*/                                          \
  X(hllc_value_type, CFL) /*Условие Куранта – Фридрихса – Леви*/          \
  X(hllc_value_type, h_min) /*шаг по пространству*/                                     \
  X(hllc_value_type, save_timer) /*интервал времени сохранения решения*/ \
  X(hllc_value_type, tau) /*шаг по времени*/

struct hllc_value_t {
#define X(type, name) type name;
  HLLC_VALUE_T
#undef X
};
static void to_json(nlohmann::json& j, const hllc_value_t& p) {
#define X(type, name) {#name, p.name},
  j = nlohmann::json{HLLC_VALUE_T};
#undef X
}
static void from_json(const nlohmann::json& j, hllc_value_t& p) {
#define X(type, name) j.at(#name).get_to(p.name);
  HLLC_VALUE_T
#undef X
}
#undef HLLC_VALUE_T

#define SOLVE_MODE_T                                                                            \
  X(int, class_vtk) /*тип сетки (и данных на ней)*/                         \
  X(int, max_number_of_iter) /*предел итераций по излучению*/          \
  X(double, accuracy) /*критерий сходимости задачи излучения*/ \
  X(int, start_point) /*начальная итерация газовой динамики*/

struct solve_mode_t {
#define X(type, name) type name;
  SOLVE_MODE_T
#undef X

  solve_mode_t() {
    class_vtk = 0;
    max_number_of_iter = 1;
    accuracy = 1e-5;
  }
};
static void to_json(nlohmann::json& j, const solve_mode_t& p) {
#define X(type, name) {#name, p.name},
  j = nlohmann::json{SOLVE_MODE_T};
#undef X
}
static void from_json(const nlohmann::json& j, solve_mode_t& p) {
#define X(type, name) j.at(#name).get_to(p.name);
  SOLVE_MODE_T
#undef X
}
#undef SOLVE_MODE_T

#ifdef DEBUG

typedef double Type;

#define ILLUM_TIMER         \
  X(Type, phys_time)        \
  X(Type, send1_wait_time)  \
  X(Type, rcv1_wait_time)   \
  X(Type, rcv2_wait_time)   \
  X(Type, rcv2_init_time)   \
  X(Type, send2_wait_time)  \
  X(Type, dir_time)         \
  X(Type, cuda_time)        \
  X(Type, cuda_wait_time_1) \
  X(Type, cuda_wait_time_2) \
  X(Type, norm_gather)

struct illum_timer_t {
#define X(type, name) type name;
  ILLUM_TIMER
#undef X
};
static void to_json(nlohmann::json& j, const illum_timer_t& p) {
#define X(type, name) {#name, p.name},
  j = nlohmann::json{ILLUM_TIMER};
#undef X
}
static void from_json(const nlohmann::json& j, illum_timer_t& p) {
#define X(type, name) j.at(#name).get_to(p.name);
  ILLUM_TIMER
#undef X
}
#undef ILLUM_TIMER

#endif  // DEBUG

#endif  // JSON_STRUCT