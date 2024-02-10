/**
 * @file global_value.h
 * @brief Файл содержит глобальные константы
 */

#ifndef GLOBAL_VALUE
#define GLOBAL_VALUE

#include <iostream>
#include <string>

#include "geo_types.h"
#include "global_consts.h"
#include "global_def.h"
#include "prj_config.h"
#include "solvers_struct.h"
#include "json/json_struct.h"

#include "log_global_consts.h"
#define LOG(val) log_##val

extern global_files_t glb_files;
extern TableFunc t_cooling_function; ///< табличная функция охлаждения

//! (в случае сплошной области задаётся так, чтобы сфера не пересекала расчётную область)
#if GEOMETRY_TYPE == Sphere
const Vector3 kCenterPoint(0, 0, 0);     ///< центр внутренней сферы на сетке
constexpr double kInternalRadius = 0.31; ///< радиус внутренней сферы (с запасом)
#elif (GEOMETRY_TYPE == TEST_ELLIPSE) || (GEOMETRY_TYPE == MAIN_ELLIPSE)
const Vector3 kCenterPoint(1, 0, 0);     ///< центр внутренней сферы на сетке
constexpr double kInternalRadius = 0.15; ///< радиус внутренней сферы (с запасом)
#else
const Vector3 kCenterPoint(10, 0, 0);    ///< центр внутренней сферы на сетке
constexpr double kInternalRadius = 0.12; ///< радиус внутренней сферы (с запасом)
#endif

//--------------------------Файлы управления---------------//
#define F_SET "settings_file.txt"
#define F_HLLC_SET "hllc_settings.txt"
#define F_INIT_HLLC "no.bin"
#define F_MPI_CONF "mpi_conf"

#define F_LOG "File_Logs.txt"

//------------------- Файлы трассировки ------------------//
#define F_GRAPH_BOUND_FACE "graph_bound_faces"
#define F_GRAPH_BODY_FACE "graph_cell_faces"
#define F_X0_LOC "LocX0"
#define F_X "X.bin"
#define F_STATE_FACE "state_face"
#define F_TRACE_VERTEX "vertex.bin"
#define F_TRACE_GRID "trace_grid.bin"

#define F_RAY_TRACE "ray_trace"
#define F_IMAGE_PLANE "image_plane"
#define F_GLOSS_CURVE "gloss_curve.txt"

//---------------Файлы трассировки сквозь внутреннюю границу------------//
#define F_DIST_TRY "dist_defining_faces"
#define F_ID_TRY "id_defining_faces"
#define F_RES "ResBound"

//----------------------Файлы для построения графа --------------------//
#define F_INTERNAL_BOUND "internal_boundary.txt"
#define F_INIT_BOUND "init_boundary.txt"
#define F_FACE_ID "face_id.txt"
#define F_GRAPH "graph"

//--------------------------Файлы геометрии-----------------------------//
#define F_GEO_FACES "geo_faces.bin"
#define F_GEO_CELLS "geo_cells.bin"

#define F_CENTERS "centers.bin"
#define F_NORMALS "normals.bin"
#define F_AREAS "areas.bin"
#define F_VOLUME "volume.bin"
#define F_NEIGHBOR "neighbors.bin"
#define F_CENTER_FACE "center_face.bin"
#define F_SURFACE "surface.bin"

//--------------------------Файлы сферы направлений-----------------------------//
#define F_DIRECTION_GRID "grid_direction.txt"
#define F_ADDITIONAL_DIRECTION_GRID "add_direction.txt"
#define F_DIRECTION_INTERPOLATION "direction_interpolation.bin" ///< номера ячеек на сфере направлений для переинтерполяции на картинную плоскость
#define F_SURFACE_SPHERE_DIRECTION "surface_grid_direction.bin" ///< ячейки на поверхности сфер направлений

//--------------------------Файлы решение-------------------------------//
#define F_SOLVE "Solve" // задается в настройках

#define F_DENSITY "density.bin"
#define F_PRESSURE "pressure.bin"
#define F_VELOCITY "velocity.bin"

#define F_ILLUM "Illum.bin"
#define F_ENERGY "energy.bin"
#define F_STREAM "stream.bin"
#define F_IMPULS "impuls.bin"
#define F_DIVSTREAM "divstream.bin"
#define F_DIVIMPULS "divimpuls.bin"

#define F_SCATTERING "int_scatter"

#if GEOMETRY_TYPE == TEST_ELLIPSE
#define F_ABSORPCOEF "alpha.bin"
#define F_RADLOOSERATE "Q.bin"
#else
#define F_ABSORPCOEF "AbsorpCoef.bin"
#define F_RADLOOSERATE "radEnLooseRate.bin"
#endif

#define F_COOLING_FUNC "rcf_H.dat"
#endif //! GLOBAL_VALUE
