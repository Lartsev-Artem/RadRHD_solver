/**
 * @file rhllc_bound_cond.h
 * @author Artem
 * @brief
 * @version 0.1
 * @date 2025-06-09
 *
 * @copyright Copyright (c) 2025
 *
 */
#if !defined RHLLC_BOUND_COND_H && defined RHLLC
#define RHLLC_BOUND_COND_H

#include "solvers_struct.h"

/*! \addtogroup rhllc Модуль расчета газовой динамики в релятивистской постановке
    @{
*/
namespace rhllc
{
    namespace bound
    {
        void OutSrcBound(const elem_t &cell, flux_all_t &flx);
        void InSrcBound(const elem_t &cell, flux_all_t &flx);
        void OuterSurfBound(const elem_t &cell, flux_all_t &flx);

    }; // namespace bound

    /**
     * @brief Функция задаёт условия контакта ячеек
     * @details Устанавливает или граничные условия или состояние соседней ячейки
     * @param[in] f грань
     * @param[in] cells ячейки сетки
     * @param[out] bound соседняя ячейка к текущей
     */
    void BoundConditions(const face_t &f, const std::vector<elem_t> &cells, flux_all_t &bound);
};

#endif // RHLLC_BOUND_COND_H