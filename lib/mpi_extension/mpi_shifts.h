
#ifndef MPI_SHIFTS
#define MPI_SHIFTS

#include <vector>

/**
 * @brief Вычисление числа компонент на узел
 *
 * @param[in] np размер кластера
 * @param[in] n общее число элементов
 * @param[out] send_count массив с количеством компонент на узел
 */
void GetSend(const int np, const int n, std::vector<int> &send_count);

/**
 * @brief Вычисление смещения компонент на узел
 *
 * @param[in] np размер кластера
 * @param[in] n общее число элементов
 * @param[out] disp массив с смещений
 */
void GetDisp(const int np, const int n, std::vector<int> &disp);

#endif // MPI_SHIFTS