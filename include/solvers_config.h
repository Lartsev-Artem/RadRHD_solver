/**
 * @file solvers_config.h
 * @author Artem
 * @brief Файл конфигурации решателей
 * @version 0.2
 * @date 2025-06-08
 *
 * @copyright Copyright (c) 2025
 *
 */

#if !defined SOLVE_CONFIG_H && defined SOLVERS
#define SOLVE_CONFIG_H

#include "prj_config.h"
#include "illum_config.h"
#include "hd_config.h"

/* ============== Available geometry forms =============== */
#define Cube 1
#define Step 2
#define Cone 3
#define Cone_JET 4
#define Sphere 5
#define Cylinder 6
#define TEST_ELLIPSE 7
#define MAIN_ELLIPSE 8
/*=======================================================*/
#define GEOMETRY_TYPE Cone

/* ============== Available task types =============== */
#define GRB_TASK 1
#define PHL_TASK 2
/*=======================================================*/
#define TASK_TYPE GRB_TASK

/* ============== Check config project =============== */

#if defined RAD_RHD && (!defined RHLLC || !defined ILLUM)
#error "Bad config. RAD_RHD is only available with RHLLC and ILLUM"
#endif
#if defined RAD_RHD && defined SPECTRUM
#error "Bad config. There can be only one task at a time (SPECTRUM or RAD_RHD)"
#endif

#endif //! SOLVE_CONFIG_H