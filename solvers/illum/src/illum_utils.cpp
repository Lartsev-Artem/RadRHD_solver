#if defined ILLUM && defined SOLVERS

#include "illum_utils.h"

#include "dbgdef.h"

/**
 * @brief Точное решение уравнения переноса (вдоль луча) в ячейке
 *
 * @param s путь
 * @param Q источник
 * @param S интеграл рассеяния
 * @param I_0 начальное значение з.Коши
 * @param k коэффициент поглощения (alpha + betta)
 * @param betta коэффициент рассеяния
 * @return значение луча в конце отрезка интегрирования
 * @note осуществляется предельный переход для среды без поглощения
 */
static inline Type GetI(Type s, Type Q, Type S, Type I_0, Type k, Type betta) {
  if (k > 1e-10)
    return (exp(-k * s) * (I_0 * k + (exp(k * s) - 1) * (Q + S * betta))) / k;
  else
    return (1 - s * k) * (I_0 + s * (Q + S * betta));
}

Type illum::BoundaryConditions(const int type_bound) {
  switch (type_bound) {
  case e_bound_free:
    return 0;

  case e_bound_lock:
    return 0;

  case e_bound_out_source:
    return 0;

  case e_bound_inner_source: {
    return 1;
#if 0
		id_try_pos++;
		grid[num_cell].nodes_value[num_in_face] = Vector3(res_on_inner_bound, res_on_inner_bound, res_on_inner_bound);
		return res_on_inner_bound;

		Type data = res_inner_bound[ShiftRes + pos_in_res++]; // защита на выход из диапазона??
		if (data >= 0) //данные от пересечения с диском или шаром
		{
			/*
				 Проверить порядок. Данные в массиве должны лежать в соответствии с упорядоченным графом
				 от направления к направлению
			*/
			return res_on_inner_bound; // I_x0;
			return data;
		}
		else // определяющимм являются противолежаащие грани (возможен расчет с учетом s=(x-x0).norm())
		{

			// результат не вполне понятен. Пока лучше использовать константу или другие параметры области (шар == граница)

			//+dist_try_surface			
			int id = id_try_surface[ShiftTry + id_try_pos - 1];  // будет лежать id грани			
			const int cell = id / 4;
			const int face = id % 4;

			Vector3 coef = grid[cell].nodes_value[face];
			Vector2	x0_local = X0[ShiftX0 + posX0++];//grid[num_cell].x0_loc[num_in_face_dir];

			I_x0 = x0_local[0] * coef[0] + x0_local[1] * coef[1] + coef[2];

			//if (I_x0 < 0) I_x0 = 0;
			return  res_on_inner_bound; // I_x0;

		}
#endif
  }
  default:
    D_LD;
  }
}

Type illum::GetIllum(const Vector3 x, const Type s, const Type I_0, const Type int_scattering, const elem_t &cell) {
  switch (_solve_mode.class_vtk) {

  case e_grid_cfg_default:
  case e_grid_cfg_static_illum: {
#if GEOMETRY_TYPE == Sphere
    Type Q = 1;
    Type alpha = 2;
    Type betta = 1;
    Type S = int_scattering;

    if ((x - Vector3(0, 0, 0)).norm() > 0.09) // излучающий шар
    {
      Q = 0;
      alpha = 0.5;
      betta = 0.5;
    }
#endif
#if GEOMETRY_TYPE == Cone
    Type Q = 0;
    Type alpha = 0.5;
    Type betta = 0.5;
    Type S = int_scattering;

    if (x[0] < 0.06) // излучающий слой
    {
      Q = 10;
      alpha = 1;
      betta = 2;
    }
#endif

    return std::max(0.0, GetI(s, Q, S, I_0, alpha + betta, betta));
  }

  case e_grid_cfg_radiation: // test task
  {
#if GEOMETRY_TYPE == TEST_ELLIPSE
    Type S = int_scattering;
    Type Q = cell.illum_val.rad_en_loose_rate; // Q=alpha*Ie
    Type alpha = cell.illum_val.absorp_coef;
    Type betta = alpha / 2; // просто из головы

    return std::max(0.0, GetI(s, Q, S, I_0, alpha + betta, betta));
#elif GEOMETRY_TYPE == MAIN_ELLIPSE

    const Type ss = s * 388189 * 1e5;          // числа --- переход к размерным параметрам
    Type Q = cell.illum_val.rad_en_loose_rate; // Q=alpha*Ie
    Type alpha = cell.phys_val.d * cell.illum_val.absorp_coef;

    return std::max(0.0, GetI(ss, Q, 0, I_0, alpha, 0));
    // I = I_0 * exp(-alpha * ss) + Q * (1 - exp(-alpha * ss)) / alpha;
    // I = I_0 * (1 - alpha * ss) + Q * ss;
#else
    D_LD;
#endif
  }

  case e_grid_cfg_full_init: // HLLC + Illum для конуса
  {
    // переход к размерным параметрам
    Type S = int_scattering * kRadiation;
    Type d = cell.phys_val.d * kDensity;
    Type p = cell.phys_val.p * kPressure;
    Type T = p / (d * kR_gas);

#if GEOMETRY_TYPE == Cone
    if (x[0] < 0.05) { //источник
      d = 0.1;
      p = 0.01;
      T = p / (d * R);
    }
#endif

    Type Ie = 2 * pow(k_boltzmann * PI * T, 4) / (15 * kH_plank * kH_plank * kH_plank * kC_Light * kC_Light);
    Type betta = kSigma_thomson * d / kM_hydrogen;
    Type alpha = betta; // alpha = ???

    Type Q = alpha * Ie;
    Type ss = s * kDist;
    Type I0 = I_0 * kRadiation;

    return std::max(0.0, GetI(ss, Q, S, I0, alpha + betta, betta) / kRadiation);
  }

  default:
    D_LD;
  }
}

Type illum::ReCalcIllum(const int num_dir, const std::vector<Vector3> &inter_coef, grid_t &grid) {

#ifdef USE_CUDA
  D_LD;
#endif

  Type norm = -1;
  const int shift_dir = num_dir * grid.size;

  for (uint32_t num_cell = 0; num_cell < grid.size; num_cell++) {

    for (int i = 0; i < CELL_SIZE; i++) {

      Vector3 Il = inter_coef[num_cell * CELL_SIZE + i];
      const Type curI = (Il[0] + Il[1] + Il[2]) / 3; //  среднее на грани (в идеале переход к ax+by+c)
      const int id = CELL_SIZE * (shift_dir + num_cell) + i;

      norm = std::max(norm, fabs(grid.Illum[id] - curI) / curI);
      grid.Illum[id] = curI;

      grid.cells[num_cell].illum_val.illum[num_dir * CELL_SIZE + i] = curI; //на каждой грани по направлениям
    }
    // cells[num_cell].illum_val.illum[num_dir] /= CELL_SIZE; //пересчёт по направлению для расчёта энергий и т.д.
  }

  return norm;
}

#endif //! defined ILLUM && defined SOLVERS