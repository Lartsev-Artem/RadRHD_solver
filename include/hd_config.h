/**
 * @file hd_config.h
 * @author Artem
 * @brief Файл конфигурации газодинамических решателей
 * @version 0.1
 * @date 2025-06-08
 *
 * @copyright Copyright (c) 2025
 *
 */
#if !defined HD_CONFIG_H && (defined HLLC || defined RHLLC)
#define HD_CONFIG_H

#include "prj_config.h"

/* ============== Available EOS =============== */
#define IDEAL 1
#define TAUB 2
/*=======================================================*/
#define EOS TAUB

/* ============== Check config project =============== */
#if defined HLLC && defined RHLLC
#error "Bad config. There can be only one task at a time (HLLC or RHLLC)"
#endif

#endif //!