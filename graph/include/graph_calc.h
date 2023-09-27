/**
 * @file graph_calc.h
 * @brief Функции расчёта фронта и формирования нового шага
 */

#if !defined GRAPH_CALC_H && defined BUILD_GRAPH
#define GRAPH_CALC_H

#include "geo_types.h"
#include "global_types.h"
#include <map>
#include <set>

/*! \addtogroup graph Модуль построения графов
    @{
*/

namespace graph {

/**
 * @brief Функция определяет текущий фронт определённой границы для сплошной сетки
 *
 * @details Проходим по всем кандидатам,
 * если число входящих сравнялось с числом определённых,
 * значит ячейка входит в текущий фронт границы.
 * @param[in] next_step_el кандидаты на вхождение в новую границу
 * @param[in] count_in_face массив с числом входящих граней ячеек
 * @param[in] count_def_face массив с числом определённых граней ячеек
 * @param[out] cur_front новая определённая граница
 * @return int ::e_type_completion
 */
int FindCurFront(const std::set<IntId> &next_step_el,
                 const std::vector<State> &count_in_face,
                 const std::vector<State> &count_def_face,
                 std::vector<IntId> &cur_front);

/**
 * @brief Функция определяет текущий фронт определённой границы для сетки с вырезанной областью
 *
 * @details Проходим по всем кандидатам,
 * если число входящих сравнялось с числом определённых,
 * значит ячейка входит в текущий фронт границы, кроме того
 * если массив outer_part не пуст, проводим пере трассировку и пытаемся определить
 * границу через геометрические объекты вне области
 * @param[in] direction направление
 * @param[in] normals нормали
 * @param[in] inner_faces грани внутренней границы
 * @param[in] next_candidates кандидаты на вхождение в новую границу
 * @param[in] count_in_face массив с числом входящих граней ячеек
 * @param[in] inner_part граничные ячейки, грани которых выходящие
 * @param[out] cur_front новая определённая граница
 * @param[out] count_def_face массив с числом определённых граней ячеек
 * @param[out] outer_part  граничные ячейки, грани которых входящие (после обнаружения, ячейка)
 * @return int ::e_type_completion
 */
int FindCurFrontWithHole(const Vector3 &direction,
                         const std::vector<Normals> &normals,
                         const std::map<IntId, FaceCell> &inner_faces,
                         const std::set<IntId> &next_candidates,
                         const std::vector<State> &count_in_face,
                         const std::set<IntId> &inner_part,
                         std::vector<IntId> &cur_front,
                         std::vector<State> &count_def_face,
                         std::set<IntId> &outer_part);

/**
 * @brief Функция формирует новый список кандидатов на определение
 *
 * @details Новая ячейка добавляется по принципу изменения соседей
 * т.е. новый потенциальный фронт формируется из всех неопределённых
 * соседей предыдущего
 * @param[in] neighbours список соседей
 * @param[in] count_in_face массив с числом входящих граней ячеек
 * @param[in] cur_front новая определённая граница
 * @param[out] count_def_face массив с числом определённых граней ячеек
 * @param[out] next_candidate кандидаты на вхождение в новую границу
 */
void NewStep(const std::vector<IntId> &neighbours,
             const std::vector<State> &count_in_face,
             const std::vector<IntId> &cur_front,
             std::vector<State> &count_def_face,
             std::set<IntId> &next_candidate);
} // namespace graph
#endif //! BUILD_GRAPH