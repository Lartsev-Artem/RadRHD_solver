#ifndef GLB_TYPES
#define GLB_TYPES

#include "global_def.h"
#include <Eigen/Dense>
#include "json/json_struct.h"

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

typedef uint8_t ShortId;
typedef const std::string& file_name_t;


enum e_boundary_types_t
{
	e_bound_free = -1,
	e_bound_lock = -2,
	e_bound_out_source = -3,
	e_bound_inner_source = -4
};

/**
 * @brief конфигурация исходной сетки
 * 
 */
enum e_grid_vtk_config_t
{
	e_grid_cfg_default 		= 0,	///< только геометрия
	e_grid_cfg_radiation 	= 1,	///< геометрия + параметры излучения излучение
	e_grid_cfg_full_init	= 2,	///< геометрия + излучение + газ
	e_grid_cfg_static_illum = 10,	///< газодинамический расчёт с однократно рассчитанным излучением
};

/**
 * @brief код завершения процесса
 * 
 */
enum e_type_completion
{
	e_completion_success = 0,	///< процесс завершился успешно
	e_completion_fail = 1		///< процесс завершился с ошибкой
};


#endif//GLB_TYPES