/**
 * @file radRHD_config.h
 * @brief Конфигурация решателя rhllc с учетом излучения
 */
#if !defined RADRHD_CONFIG_H && defined RAD_RHD
#define RADRHD_CONFIG_H

#define ISOTROPIC_INTENSITY 0 ///< подключает изотропный расчёт давления излучения через тензор эддингтона (не меняет расчёт illum модуля)

#endif //! RADRHD_CONFIG_H