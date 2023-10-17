
#if !defined HLLC_FLUX_H && defined SOLVERS && defined HLLC
#define HLLC_FLUX_H

/*! \addtogroup hllc Модуль расчета газовой динамики
    @{
*/

namespace hllc {

/**
 * @brief Расчёт консервативных переменных через физические
 *
 * @details статья ost1098
 * @param[in] W физические переменные
 * @param[out] U консервативные переменные
 */
void GetConvValue(const flux_t &W, flux_t &U);

/**
 * @brief Восстановление физических переменных через консервативные
 *
 * @details алгоритм построен на основе метода Ньютона. из статьи ost1098
 * @param[in] U консервативные переменные
 * @param[out] W физические переменные
 */
void GetPhysValue(const flux_t &U, flux_t &W);

/**
 * @brief Расчёт потоков методом HLLC
 *
 * @param[in] conv_val_l консервативный ячейка слева
 * @param[in] conv_val_r консервативный ячейка справа
 * @param[in] phys_val_l физическая ячейка слева
 * @param[in] phys_val_r физическая ячейка справа
 * @param[out] f поток
 */
void GetFlux(const flux_t &conv_val_l, const flux_t &conv_val_r,
             const flux_t &phys_val_l, const flux_t &phys_val_r, face_t &f);

} // namespace hllc

#endif //! HLLC_FLUX_H