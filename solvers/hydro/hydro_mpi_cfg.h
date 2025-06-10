#if !defined HYDRO_MPI_CFG_H && USE_MPI && (defined RHLLC || defined HLLC)
#define HYDRO_MPI_CFG_H

#include "solvers_struct.h"

namespace rrhd
{
    namespace hydro_mpi
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

    } // namespace hydro_mpi
} // namespace rrhd

#endif //! HYDRO_MPI_CFG_H