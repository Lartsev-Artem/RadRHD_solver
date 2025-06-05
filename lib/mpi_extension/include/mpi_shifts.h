/**
 * @file mpi_shifts.h
 * @brief Расчет локальных смещений в линейных массивах на разных узлах
 *
 */
#ifndef MPI_SHIFTS
#define MPI_SHIFTS

#include "global_types.h"
#include <vector>

#include "mpi_ext.h"
/**
 * @brief Вычисление число компонент на узел
 *
 * @param[in] np размер кластера
 * @param[in] n общее число элементов
 * @param[out] send_count массив с количеством компонент на узел
 */
void GetSend(const int np, const IdType n, std::vector<IdType> &send_count);

/**
 * @brief Вычисление смещения компонент на узел
 *
 * @param[in] np размер кластера
 * @param[in] n общее число элементов
 * @param[out] disp массив с смещений
 */
void GetDisp(const int np, const IdType n, std::vector<IdType> &disp);

#ifdef USE_MPI
/**
 * @brief Устанавливает число компонент на узел и смещения компонент на узел
 *
 * @param[in] metis_id карта ячеек по узлам (важно кол-во вхождений на узел а не порядок(т.е. сетка может быть перенумерованной))
 * @param[in] np число процессов
 * @param[inout] send размер кластеров
 * @param[inout] disp смешение на каждам из узлов
 */
void SetShifts(const std::vector<int> &metis_id, int np,
               std::vector<int> &send, std::vector<int> &disp);

#endif
#endif // MPI_SHIFTS