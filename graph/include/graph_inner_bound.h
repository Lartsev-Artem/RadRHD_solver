/**
 * @file graph_inner_bound.h
 * @brief Поиском пересечений на внутренний границе
 *
 * @warning требует поддержку cuda
 */
#if !defined GRAPH_INNER_BOUND_H && defined BUILD_GRAPH
#define GRAPH_INNER_BOUND_H

#include "geo_types.h"
#include <map>
#include <set>

/*! \addtogroup graph Модуль построения графов
    @{
*/

namespace graph {

/**
 * @brief Подпространство имён построителя графов с поиском пересечений на внутренний границе
 *
 */
namespace trace_through_boundary {

struct IntersectBound_t {
  std::vector<IntId> code;        ///<  коды пересечений
  std::vector<IntId> out_id_cell; ///< номера ячеек на выходящей границе
};

/**
 * @brief Инициализирует память на видеокарте и копирует геометрию границы
 * @note память выделяется под всю границу целиком
 *
 * @return int ::::e_type_completion
 */
int InitDevice();

void ClearDevice();

/**
 * @brief Поиск ячеек и пересечений с геометрией на внутренней границе
 *
 * @details Функция находит коды пересечений на внутренней границе
 * код ячейки --- глобальный номер грани, код объекта ::e_ray_intersect_code
 * @param[in] num_dir номер направления в сфере направлений
 * @param[in] direction прямое направление трассировки
 * @param[in] faces полная граница с ключом - номером ячейки
 * @param[in] outer_part номера ячеек на выходящей границе (по ним формируется матрица трассировки лучей)
 * @param[out] intersections коды пересечений (не упорядоченные с графов)
 * @return int ::e_type_completion
 * @warning видеокарта должна быть проинициализирована
 */
int FindBoundCondOnInnerBoundary(int num_dir, const Vector3 &direction, const std::map<IntId, FaceCell> &faces, const std::set<IntId> &outer_part,
                                 std::vector<IntId> &intersections);

/**
 * @brief Переупорядочивание ячеек в соответствие с графом направления для конвеерного чтения в решателе
 *
 * @param[in] graph граф направление
 * @param[inout] intersections структура связывающая коды и номера ячеек
 */
void SortInnerBoundary(const std::vector<IntId> &graph, IntersectBound_t &intersections);

} // namespace trace_through_boundary
} // namespace graph
#endif //! GRAPH_INNER_BOUND_H