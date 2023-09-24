/**
 * @file graph_config.h
 * @brief Файл конфигурации построителя графов
 *
 */

#if !defined GRAPH_CONFIG_H && defined BUILD_GRAPH
#define GRAPH_CONFIG_H

#include "prj_config.h"

//#define USE_OMP        ///< подключение технологии omp (самостоятельно / вместе mpi)

#define GRID_WITH_INNER_BOUNDARY ///< граф для сетки с внутренней границей (для оптимизации следует отключить, при использование сплошной сетки)

//#define USE_STRANGE_FUNCTION ///< какие то старые функции, которые в текущем конфиге не используются

#ifdef DEBUG
//#define ONLY_ONE_DIRECTION  // только в последовательном режиме
#endif

#endif //! GRAPH_CONFIG_H