/**
 * @file graph_config.h
 * @brief Файл конфигурации построителя графов
 *
 */

#if !defined GRAPH_CONFIG_H && defined BUILD_GRAPH
#define GRAPH_CONFIG_H

/*! \addtogroup graph Модуль построения графов
    @{
*/

#include "prj_config.h"

// #define GRAPH_TRACING_INNER_BOUNDARY ///< включать трассировку сквозь внутреннюю границу
#define GRID_WITH_INNER_BOUNDARY ///< граф для сетки с внутренней границей (для оптимизации следует отключить, при использование сплошной сетки)

#ifdef RRHD_DEBUG
// #define ONLY_ONE_DIRECTION  // только в последовательном режиме
#endif

#endif //! GRAPH_CONFIG_H