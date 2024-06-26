/**
 * @file graph_init_state.h
 * @brief Функции инициализации структур по выбранному направлению
 * @bug для сеток с внутренний границей алгоритм может завершаться
 * до окончания построения графа. В качестве исправления предлагается вызывать функцию
 * ::TryRestart, которая заново инициализирует текущую границу и запускает алгоритм с
 * точки завершения. число рестартов ограничено
 */

#if defined BUILD_GRAPH
#ifndef GRAPH_INIT_STATE_H
#define GRAPH_INIT_STATE_H

#include "geo_types.h"
#include "global_types.h"

#include <map>
#include <set>

/*! \addtogroup graph Модуль построения графов
    @{
*/

namespace graph {

constexpr int GRAPH_MAX_RESTART = 5; ///< максимальное число рестартов на направление

/**
 * @brief Функция запускает процесс рестарта поиска текущей границы в случае преждевременной остановки
 *
 * @param[in] count_in_face массив количества входящих граней ячейки
 * @param[in] count_def_face массив количества определённых граней ячейки
 * @param[in] outer_part номера ЯЧЕЕК выходящей части
 * @param[out] cur_front новая определённая граница
 * @param[out] next_candidate начальный набор ячеек-кандидатов для упорядоченного графа
 * @return int ::e_type_completion
 */
int TryRestart(const std::vector<State> &count_in_face,
               const std::vector<State> &count_def_face,
               const std::set<IntId> &outer_part,
               std::vector<IntId> &cur_front,
               std::set<IntId> &next_candidate);

/**
 * @brief Функция определяет состояние граней по отношению к границе
 *
 * \details состояние undef- означает, что грань требует определения через другую.
 * def - грань может быть определена (на момент инициализации через гран. условия)
 * \note внутренняя граница требует определения. её состояние задаётся undef
 * @param[in] neighbours все соседи граней
 * @param[in] inter_faces грани внутренней границы
 * @param[out] faces_state массив состояний ::e_face_state_t
 */
void InitFacesState(const std::vector<IntId> &neighbours,
                    const std::map<IntId, FaceCell> &inter_faces,
                    std::vector<State> &faces_state);

/**
 * @brief Функция разделяет граничные ячейки на входящие и выходящие в зависимости от направления
 *
 * @note выходящая ячейка требует трассировки через внутреннюю область
 * выходящая --- значит луч входит в нее из внутренний области
 * входящая --- значит луч через нее выходит из основной области
 * @param[in] direction выбранное направление
 * @param[in] normals нормали к ячейкам
 * @param[in] inter_boundary_face_id глобальная нумерация  граней внутренней границы
 * @param[out] inner_part номера ЯЧЕЕК входящей части
 * @param[out] outer_part номера ЯЧЕЕК выходящей части
 */
void DivideInnerBoundary(const Vector3 &direction,
                         const std::vector<Normals> &normals,
                         const std::set<IntId> &inter_boundary_face_id,
                         std::set<IntId> &inner_part,
                         std::set<IntId> &outer_part);

// число входящих граней для каждой ячейки + число известных из них + начальная граница

/**
 * @brief Определение числа входящих и определённых граней, формирующих начальный набор графа
 *
 * @param[in] dir направление
 * @param[in] normals нормали к ячейкам
 * @param[in] faces_state массив состояний ::e_face_state_t
 * @param[out] count_in_face массив количества входящих граней ячейки
 * @param[out] count_def_face массив количества определённых граней ячейки
 * @param[out] next_step_el начальный набор ячеек-кандидатов для упорядоченного графа
 */
void FindNumberOfAllInAndDefFaces(const Vector3 &dir,
                                  const std::vector<Normals> &normals,
                                  const std::vector<State> &faces_state,
                                  std::vector<State> &count_in_face, std::vector<State> &count_def_face,
                                  std::set<IntId> &next_step_el);

} // namespace graph
#endif //! GRAPH_INIT_STATE_H
#endif // BUILD_GRAPH