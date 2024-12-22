/**
 * @file ray_tracing_const.h
 * @brief Константы с описанием геометрии
 *
 */
#ifndef RAY_TRACING_CONST_H
#define RAY_TRACING_CONST_H

#include "geo_types.h"

/*! \addtogroup ray_tracing Модуль трассировки лучей, построения картинной плоскости и кривых блеска
    @{
*/

/**
 * @brief Код объекта пересечения
 * \note неотрицательные значения будут соответствовать ячейкам (или граням)
 *
 */
enum e_ray_intersect_code {
  e_ray_intersect_none = -1,   ///< нет пересечения
  e_ray_intersect_disk = -2,   ///< пересечение с диском
  e_ray_intersect_sphere = -3, ///< пересечение со сферой
  e_ray_intersect_rosh = -4    ///< пересечение с потенциалом роша
};

#ifdef __NVCC__
#define HOST_DEVICE __device__
#else
#define HOST_DEVICE
#endif

namespace ray_tracing {

HOST_DEVICE constexpr double PI = 3.1415926535897932384626433832795; ///<число пи

HOST_DEVICE constexpr Type k_center_x = 1, k_center_y = 0, k_center_z = 0; ///< координаты центр аккретора
const Vector3 k_center_sphere(k_center_x, k_center_y, k_center_z);         ///< центр аккретора
HOST_DEVICE constexpr Type k_radius_sphere = 0.01;                         // 0.01;                         ///< радиус аккретора
HOST_DEVICE constexpr Type k_internal_radius_disk = 0.01;                  ///< внутренний радиус около аккреционного диска
HOST_DEVICE constexpr Type k_external_radius_disk = 0.06;                  ///< внешний радиус около аккреционного диска

// параметры картинной плоскости
#ifdef GRB_TASK
constexpr Type k_width_plane = 1.2;  ///< безразмерная ширина картинной плоскости
constexpr Type k_height_plane = 0.6; ///< безразмерная высота картинной плоскости
#else
constexpr Type k_width_plane = 2.4;  ///< безразмерная ширина картинной плоскости
constexpr Type k_height_plane = 1.2; ///< безразмерная высота картинной плоскости
#endif

constexpr int k_pixels_width = 600;         ///< число пикселей в ширину
constexpr int k_pixels_height = 300;        ///< число пикселей в высоту
constexpr Type k_height_above_center = 0.1; ///< высота над центром масс

const Vector3 k_center_of_mass(0.879518, 0, 0); ///< центр масс системы
constexpr int k_number_of_frame = 120;          ///< число позиций в окружности

constexpr Type k_accretion_energy = 10; ///< энергия на поверхности аккретора
constexpr Type k_disk_energy = 5;       ///< энергия аккреционного диска около аккретора
constexpr Type k_rosh_energy = 0.00001; ///< энергия на поверхности потенциала Роша

/**
 * @brief Параметры картинной плоскости в 2d координатах
 */
struct PlaneParams
{
public:
  Type _width; ///< безразмерная ширина плоскости
  Type _height;///< безразмерная высота плоскости
  int _pixels_width; ///< число пикселей в ширину
  int _pixels_height; ///< число пикселей в высоту
private:
  Vector3 _angle_of_plane; ///< угол плоскости. От него начинается заполнение всей плоскости
  Type _step_x;  ///< ширина пикселя
  Type _step_y;  ///< высота пикселя
public:
  PlaneParams() = delete;
  PlaneParams(const PlaneParams& prm)
  {
    *this = PlaneParams(prm._width, prm._height, prm._pixels_width, prm._pixels_height);
  }
  PlaneParams(Type width,Type height,int pixels_width, int pixels_height) :
  _width(width),
  _height(height),
  _pixels_width(pixels_width),
  _pixels_height(pixels_height) 
  {
    _angle_of_plane = Vector3(-(_width / 2), -(_height / 2), 0);
    _step_x = _width / _pixels_width;
    _step_y = _height / _pixels_height;
  }
  void operator=(const PlaneParams& prm)
  {
    _width = prm._width;
   _height = prm._height;
   _pixels_width=prm._pixels_width;
   _pixels_height = prm._pixels_height;

    _angle_of_plane = prm._angle_of_plane;
    _step_x=prm._step_x;
    _step_y=prm._step_y;
  }

  /**
   * @brief Возвращает координаты центра [i,j] пикселя
   * 
   * @param[in] i номер пикселя по горизонтали
   * @param[in] j номер пикселя по вертикали
   * @return Vector3 
   */
  Vector3 get_pixel_coord(const int i, const int j) const
  {
      return Vector3(_angle_of_plane(0) + i * _step_x, _angle_of_plane(1) + j * _step_y, 0);
  }
};

/**
 * @brief Параметры сцены для проецирования
 */
struct ParamTraceProjection
{
  ParamTraceProjection() = delete;
  ParamTraceProjection(const PlaneParams& prm2D, const Vector3& orig, const Vector3& observe_dir)
                      :params2D(prm2D),plane_orig(orig),observer(observe_dir) {}
  
  PlaneParams params2D; ///< параметры плоскости
  Vector3 plane_orig; ///< центр плоскости
  
  //для разных моделей направление лучей может быть параллельно одному
  // направлению или проходить через одну точку
  union 
  {
    Vector3 observer; ///< положение наблюдателя
    Vector3 direction; ///< направление лучей
  };
};

} // namespace ray_tracing
#endif //! RAY_TRACING_CONST_H