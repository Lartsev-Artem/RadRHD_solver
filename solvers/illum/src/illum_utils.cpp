#if defined ILLUM && defined SOLVERS

#include "illum_utils.h"

#include "dbgdef.h"

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

#endif //! defined ILLUM && defined SOLVERS