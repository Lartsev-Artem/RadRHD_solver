/**
 * @file global_types.h
 * @brief Файл с глобальными определениями (и псевдонимами) типов данных
 *
 */
#ifndef GLB_TYPES
#define GLB_TYPES

#include <Eigen/Dense>

//#include "global_def.h"
//#include "json/json_struct.h"

typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector2d Vector2;
typedef Eigen::VectorXd VectorX;
typedef Eigen::Matrix3d Matrix3;
typedef Eigen::Vector4d Vector4;
typedef Eigen::Matrix4d Matrix4;
typedef Eigen::MatrixXd MatrixX;

typedef double Type;
typedef int IntId;
typedef uint8_t State;
typedef uint32_t bits_flag_t; ///< битовые флаги

typedef uint8_t ShortId;
typedef const std::string &file_name_t;

/**
 * @brief тип граничные условий
 *
 */
enum e_boundary_types_t {
  e_bound_free = -1,        ///< свободная граница
  e_bound_lock = -2,        ///< стенка
  e_bound_out_source = -3,  ///< внешний источник
  e_bound_inner_source = -4 ///< внутренний источник
};

/**
 * @brief конфигурация исходной сетки
 *
 */
enum e_grid_vtk_config_t {
  e_grid_cfg_default = 0,       ///< только геометрия
  e_grid_cfg_radiation = 1,     ///< геометрия + параметры излучения излучение
  e_grid_cfg_full_init = 2,     ///< геометрия + излучение + газ
  e_grid_cfg_static_illum = 10, ///< газодинамический расчёт с однократно рассчитанным излучением
};

/**
 * @brief код завершения процесса
 *
 */
enum e_type_completion {
  e_completion_success = 0, ///< процесс завершился успешно
  e_completion_fail = 1     ///< процесс завершился с ошибкой
};

/**
 * @brief код cuda потока
 *
 */
enum e_cuda_stream_id_t {
  e_сuda_scattering_1 = 0, ///< первая часть интеграла рассеяния
  e_сuda_scattering_2 = 1, ///< вторая часть интеграла рассеяния
  e_сuda_params = 2,       ///<  расчёт физических величин (потоки, импульсы)
  // e_сuda_sender,
  e_сuda_count ///< общее число потоков
};

#endif // GLB_TYPES