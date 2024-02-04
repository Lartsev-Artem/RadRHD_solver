/**
 * @file compton.h
 * @brief Функции расчёта комптон эффектов
 * @version 0.1
 * @date 2024-02-04
 * @warning Скорости размерные!!!
 *
 */
#ifndef COMPTON_H
#define COMPTON_H

#include "geo_types.h"

/**
 * @brief Функция вычисляет полное сечение рассеяния
 *
 * @param frq частота фотона
 * @return полное сечение рассеяния
 */
double get_scat_coef(double frq);

/**
 * @brief Функция вычисляет полное сечение рассеяния в релятивистском приближении
 *
 * @param frq частота
 * @param vel модуль скорости
 * @param cosf косинус между направлениями электрона и фотона до взаимодействия
 * @return полное сечение рассеяния
 */
double get_scat_coef(double frq, double vel, double cosf);

/**
 * @brief Функция вычисляет комптоновскую частоту в релятивистском приближении
 *
 * @param frq частота
 * @param vel модуль скорости
 * @param cosf косинус между направлениями электрона и фотона до взаимодействия
 * @param cosf1 косинус между направлениями электрона и фотона после взаимодействия
 * @param cosTh косинус между направлениями фотона  до взаимодействия и после
 * @return комптоновская частота
 */
double get_compton_frq(double frq, double vel, double cosf, double cosf1, double cosTh);

/**
 * @brief Функция вычисляет дифференциальное сечение рассеяния в релятивистском приближении
 *
 * @param frq частота
 * @param vel модуль скорости
 * @param cosf косинус между направлениями электрона и фотона до взаимодействия
 * @param cosf1 косинус между направлениями электрона и фотона после взаимодействия
 * @param cosTh косинус между направлениями фотона  до взаимодействия и после
 * @return double
 */
double get_dif_scat_coef(double frq, double vel, double cosf, double cosf1, double cosTh);

/**
 * @brief Функция вычисляет подинтегральную функцию рассеяния
 *
 * @param frq частота
 * @param vel скорость
 * @param dir направление до взаимодействия
 * @param dir_scat направление после взаимодействия
 * @param Illum массив излучения по дискретным частотам в данной точке и направлении
 * @return значение под интегралом рассеяния
 */
Type get_int_func(const Type frq, const Vector3 &vel, const Vector3 &dir, const Vector3 &dir_scat, const Type *Illum);

#endif //! COMPTON_H