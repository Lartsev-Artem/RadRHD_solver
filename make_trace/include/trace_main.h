/**
 * @file trace_main.h
 * @brief Файл подключает модуль трассировки для метода коротких характеристик
 * @version 0.1
 * @date 2023-09-28
 *
 * @copyright Copyright (c) 2023
 *
 */
#if !defined TRACE_MAIN_H && defined MAKE_TRACE
#define TRACE_MAIN_H

/*! \addtogroup trace Модуль предварительной трассировки
    \brief Модуляь реализует трассировку сквозь сетку по выбранным направлениям.
    \note требует построенные графы обхода сетки (модуль ::graph)
    @{
*/

/**
 * @brief Пространство имён модуля трассировки
 *
 */
namespace trace {

/**
 * @brief Функция вызывает модуль трассировки по всем направлениям
 * \note  Перед вызовом должна быть проинициализирована глобальная структура ::global_files_t
 * построена внутренняя геометрия с помощью ::BuildDataFromVTK, и построены графы обхода с помощью ::RunGraphModule
 *
 * @return int ::e_type_completion
 */
int RunTracesModule();

} // namespace trace

#endif //! TRACE_MAIN_H