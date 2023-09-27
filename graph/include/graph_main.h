/**
 * @file graph_main.h
 * @brief Файл подключает модуль построителя графов
 *
 * @details
 * @version 0.1
 * @date 2023-09-25
 *
 */
#if !defined GRAPH_MAIN_H && defined BUILD_GRAPH
#define GRAPH_MAIN_H

/*! \addtogroup graph Модуль построения графов
    \brief Модуля реализует поиск порядка обхода ячеек для маршевого алгоритма
    расчёта излучения в области
    @{
*/

/**
 * @brief Пространство имён модуля построителя упорядоченного графа
 *
 */
namespace graph {

/**
 * @brief Функция вызывает модуль построителя по всем направления графов
 * \note  Перед вызовом должна быть проинициализирована глобальная структура ::global_files_t
 * и построена внутренняя геометрия с помощью ::BuildDataFromVTK
 *
 * @return int ::e_type_completion
 */
int RunGraphModule();
} // namespace graph

#endif //! GRAPH_MAIN_H