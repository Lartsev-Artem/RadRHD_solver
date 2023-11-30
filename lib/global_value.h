/**
 * @file global_value.h
 * @brief Файл содержит глобальные константы
 * \todo физические константы
 */

#ifndef GLOBAL_VALUE
#define GLOBAL_VALUE

#include <iostream>
#include <string>

#include "geo_types.h"
#include "global_def.h"
#include "prj_config.h"
#include "solvers_struct.h"
#include "json/json_struct.h"

extern global_files_t glb_files;
extern TableFunc t_cooling_function; ///< табличная функция охлаждения

constexpr double PI = 3.1415926535897932384626433832795; ///<число пи

constexpr double kMinPressure = 1e-12;
constexpr double kMinDensity = 1e-12;

constexpr double kGamma1 = 4. / 3; ///< показатель адиабаты
constexpr double kGamma_g = kGamma1 / (kGamma1 - 1.0);

//! (в случае сплошной области задаётся так, чтобы сфера не пересекала расчётную область)
#if GEOMETRY_TYPE == Sphere
const Vector3 kCenterPoint(10, 0, 0);    ///< центр внутренней сферы на сетке
constexpr double kInternalRadius = 0.31; ///< радиус внутренней сферы (с запасом)
#elif (GEOMETRY_TYPE == TEST_ELLIPSE) || (GEOMETRY_TYPE == MAIN_ELLIPSE)
const Vector3 kCenterPoint(1, 0, 0);     ///< центр внутренней сферы на сетке
constexpr double kInternalRadius = 0.15; ///< радиус внутренней сферы (с запасом)
#else
const Vector3 kCenterPoint(10, 0, 0);    ///< центр внутренней сферы на сетке
constexpr double kInternalRadius = 0.12; ///< радиус внутренней сферы (с запасом)
#endif

#if 1 //СГС

constexpr double kC_Light = 3 * 1e10;                     ///<скорость света в м/c
constexpr double kC_LightInv = (1.0 / (kC_Light));        ///< обратная величина к скорости света
constexpr double kR_gas = 83144626.1815324;               ///< газовая постоянная [ Дж/(моль*К)]
constexpr double kH_plank = 6.62 * 1e-27;                 ///< постоянная Планка[кг * м^2 /с]
constexpr double k_boltzmann = 1.3807 * 1e-16;            ///< постоянная Больцмана[Дж/K] = [ кг*м^2/(с^2*T)]
constexpr double kSigma_thomson = 6.65210 * 1e-25;        ///< сечение томсоновского рассеяния [m^2]
constexpr double kM_hydrogen = 1.6735575 * 1e-24;         ///< масса водорода[кг]
constexpr double kDistAccretor = 3.88190065213158 * 1e10; ///< характерное расстояние
constexpr double kStefanBoltzmann = 5.670374419 * 1e-5;   ///< постоянная Стефана-Больцмана[ эрг·с^−1·см^−2·К^−4]

#else //СИ

constexpr double kC_Light = 299792458.0;                ///<скорость света в м/c
constexpr double kC_LightInv = (1.0 / (kC_Light));      ///< обратная величина к скорости света
constexpr double kR_gas = 8.314;                        ///< газовая постоянная [ Дж/(моль*К)]
constexpr double kH_plank = 6.62 * 1e-34;               ///< постоянная Планка[кг * м^2 /с]
constexpr double k_boltzmann = 1.38 * 1e-23;            ///< постоянная Больцмана[Дж/K] = [ кг*м^2/(с^2*T)]
constexpr double kSigma_thomson = 6.65210 * 1e-29;      ///< сечение томсоновского рассеяния [m^2]
constexpr double kM_hydrogen = 1.6735575 * 1e-27;       ///< масса водорода[кг]
constexpr double kStefanBoltzmann = 5.670374419 * 1e-8; ///< постоянная Стефана-Больцмана[Вт*м^-2*К^−4]
#endif

constexpr double kEarthMass = (5.9722 * 1e25); ///< масса Земли в кг
constexpr double kSunMass = (1.9891 * 1e31);   ///< масса Солнца в кг

constexpr double kDistSun = (149.6 * 10e9); ///< расстояние до Солнца в м
constexpr double kDistMoon = 400000000.;    ///<расстояние до Луны в м

constexpr double kDist = 1e6;          ///< характерное расстояние
constexpr double kMass = (1 * 1e21);   ///< характерная масса
constexpr double kVelocity = kC_Light; ///< характерная скорость
// constexpr double kTemperature = 4000;         ///< характерная температура
constexpr double kTime = (kDist / kVelocity); ///< характерное время

// constexpr double kDensity = (kMass / (kDist * kDist * kDist));   ///< характерная плотность
// constexpr double kPressure = (kMass / (kDist * kTime * kTime));  ///< характерное давление
// constexpr double kRadiation = (kMass / (kTime * kTime * kTime)); ///< характерное излучение

constexpr double kDensity = (3.34 * 10e-14);                                  ///< характерная плотность
constexpr double kPressure = (kDensity * kVelocity * kVelocity);              ///< характерное давление
constexpr double kRadiation = (kDensity * kVelocity * kVelocity * kVelocity); ///< характерное излучение

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
