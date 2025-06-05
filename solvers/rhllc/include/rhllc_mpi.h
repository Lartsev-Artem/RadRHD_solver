/**
 * @file rhllc_mpi.h
 * @brief Функции организующие обмен данных между узлами
 * @version 0.1
 * @date 2024-02-17
 *
 */
#if !defined RHLLC_MPI_H && defined SOLVERS && defined USE_MPI
#define RHLLC_MPI_H

#include "mpi_ext.h"
#include "solvers_struct.h"

/*! \addtogroup rhllc Модуль расчета газовой динамики в релятивистской постановке
    @{
*/

namespace rhllc_mpi
{

    /**
     * @brief Инициализация рассылки физических переменных всем узлам
     *
     * @param[inout] hllc_st структура mpi пересылок(должна быть инициализирована)
     * @param[inout] grid сетка
     */
    void StartPhysCast(mpi_hd_t *hllc_st, grid_t &grid);

    /**
     * @brief Синхронизация процессов (после рассылки)
     * @param[inout] hllc_st структура mpi пересылок
     */
    void SyncPhysCast(mpi_hd_t *hllc_st);

    /**
     * @brief Синхронизация процессов (после рассылки) с одновременной расчётом физических параметров для излучения
     * @param[inout] grid сетка
     */
    void SyncAndCalcPhysCast(grid_t &grid);

    /**
     * @brief Инициация рассылки (приёма передачи) физических и консервативных переменных
     * на границе узлов
     * @param[inout] hllc_st структура mpi пересылок
     */
    void StartExchangeBoundaryCells(mpi_hd_t *hllc_st);

    /**
     * @brief Синхронизация (приёма передачи) физических и консервативных переменных
     * @param[inout] hllc_st структура mpi пересылок
     */
    void SyncExchangeBoundaryCells(mpi_hd_t *hllc_st);

    /**
     * @brief Инициализации конфигурации MPI пересылок
     * @details Функция полностью инициализирует структуру mpi_hd_t.
     * Определяется геометрия и расположение узлов, инициализируются mpi запросы
     * Ожидается, что ячейки заранее сгруппированы.
     * Карта кластера имеет вид: [(c0,c1,c2), {c3,|c4}, (c5,c6,c7) ],
     * где ci - ячейки
     * |  - разделение по узлам
     * {} - ячейки подлежащие обмена
     * () - регулярные ячейки определенные на узле
     * [] - граница глобальной сетке
     *
     * @param[in] metis_id распределение ячеек по узлам
     * @param[inout] grid сетка(инициализирует данные о принадлежности к узлу)
     * @param[in] comm коммуникатор на котором производится расчет газовой динамики
     * @note должны быть инициализированы типы mpi
     * @warning Каждый узел должен иметь не более двух соседей (по пространству)

     */
    void InitMpiConfig(const std::vector<int> &metis_id, grid_t &grid, MPI_Comm comm = MPI_COMM_WORLD);

    void Hllc3dStab(const Type tau, grid_t &grid);

    void HllcConvToPhys(grid_t &grid);

    /**
     * @brief Добавление влияния излучения в газовую динамику
     * @param[inout] grid сетка с вычисленным излучением
     */
    void AddRadFlux(grid_t &grid);
} // namespace rhllc_mpi

#endif //! RHLLC_MPI_H