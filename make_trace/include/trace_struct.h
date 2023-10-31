/**
 * @file trace_struct.h
 * @brief Локальные структуры для модуля трассировки
 *
 */

#if !defined TRACE_STRUCT_H && defined MAKE_TRACE
#define TRACE_STRUCT_H

#include "global_types.h"

namespace trace {

// параметры диска и внутренней сферы:
constexpr double kRsphere = 0.001; ///< радиус аккретора
constexpr double kR1disk = 0.001;  ///< внутренний радиус диска
constexpr double kR2disk = 0.09;   ///< внешний радиус диска
} // namespace trace
#endif //! TRACE_STRUCT_H