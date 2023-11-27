/**
 * @file rhllc_flux_stab.h
 * @brief Расчёт состояние RHD. За основу взят код PLUTO.
 *
 */
#if !defined RHLLC_FLUX_STAB_H && defined SOLVERS
#define RHLLC_FLUX_STAB

#include "solvers_struct.h"

/*! \addtogroup rhllc Модуль расчета газовой динамики в релятивистской постановке
    @{
*/

namespace rhllc {

/*!
 *
 * \param [in,out]  u      array of conservative variables (entropy will
 *                         be redefined)
 * \param [out]     v      array of primitive variables
 *
 * \return Error codes:
 *  - 0 = success
 *  - 1 = solution does not exist
 *  - 2 = negative pressure
 *  - 3 = inaccurate solution (debug only)
 *  - 4 = NaN
 *********************************************************************** */

/**
 * @brief Расчёт консервативных переменных через физические с контролем состояния
 *
 * @details Осуществляется контроль максимальной скорости. В случае превышения
 * устанавливается верхний предел с пересчетом физического состояния
 * @param[inout] W физические переменные
 * @param[out] U консервативные переменные
 */
void GetConvValueStab(flux_t &W, flux_t &U);

/**
 * @brief Восстановление физических переменных через заданное минимальное давление
 *
 * @details Фиксируем давление на \c ::kMinPressure  и находим четыре-скорости,
    решая соотношение импульс-скорость (в квадрате):
    \f[
        m/(Dh) - u = 0
    \f]
    используя метод секущих, в предположении:
    - u(0) = m/D : это максимально допустимая скорость, для которой p = 0;
    - u(1) = u+ : положительный корень квадратного уравнения для которого \f$ sqrt(1+u*u) = u\f$.

 * @note Это уравнение всегда должно иметь решение.
 * @param[inout] U консервативные переменные
 * @param[out] W физические переменные
 * \return Error codes are:
 * - 0 = success
 * - 1 = v^2 > 1
 * - 2 = too many iterations
 */
int PhysPressureFix(flux_t &U, flux_t &W);

/**
 * @brief Восстановление физических переменных через консервативные
 *
 * @details Попытка восстановить физические переменные, вычислив
            давление газа с использованием общей плотности энергии,
            используя поиск корней Ньютона-Рафсона.

 * @param[in] U консервативные переменные
 * @param[out] W физические переменные
 * @return Error codes:
 *  - 0 = success
 *  - 1 = solution does not exist
 *  - 2 = negative pressure
 *  - 3 = inaccurate solution (debug only)
 *  - 4 = NaN
 */
int GetPhysValueStab(const flux_t &U, flux_t &W);

/**
 * @brief Расчёт потоков методом HLLC
 *
 * @param[in] conv_val_l консервативный ячейка слева
 * @param[in] conv_val_r консервативный ячейка справа
 * @param[in] phys_val_l физическая ячейка слева
 * @param[in] phys_val_r физическая ячейка справа
 * @param[out] f поток
 * @return максимальная скорость сигнала
 */
Type GetFluxStab(const flux_t &conv_val_l, const flux_t &conv_val_r,
                 const flux_t &phys_val_l, const flux_t &phys_val_r, face_t &f);

} // namespace rhllc

#endif //! RHLLC_FLUX_H