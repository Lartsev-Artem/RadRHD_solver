/**
 * @file global_consts.h
 * @brief Файл содержит глобальные определения всех констант
 *
 */
#ifndef GLOBAL_CONSTS_H
#define GLOBAL_CONSTS_H

#include "log_global_consts.h"
#define LOG(val) log_##val

constexpr double k_exp = 2.7182818284590452353602874713527; ///<число e
constexpr double PI = 3.1415926535897932384626433832795;    ///<число пи
constexpr double PI4 = 4 * PI;
constexpr double kMinPressure = 1e-12;
constexpr double kMinDensity = 1e-12;

constexpr double kGamma1 = 4.0 / 3; ///< показатель адиабаты
constexpr double kGamma_g = kGamma1 / (kGamma1 - 1.0);

#if 1 //СГС

constexpr double kC_Light = 3 * 1e10;                     ///<скорость света в м/c
constexpr double kC_LightInv = (1.0 / (kC_Light));        ///< обратная величина к скорости света
constexpr double kR_gas = 83144626.1815324;               ///< газовая постоянная [ Дж/(моль*К)]
constexpr double kH_plank = 6.62 * 1e-27;                 ///< постоянная Планка[кг * м^2 /с]
constexpr double k_boltzmann = 1.3807 * 1e-16;            ///< постоянная Больцмана[Дж/K] = [ кг*м^2/(с^2*T)]
constexpr double kSigma_thomson = 6.65210 * 1e-25;        ///< сечение томсоновского рассеяния [m^2]
constexpr double kM_hydrogen = 1.6735575 * 1e-24;         ///< масса водорода[кг]
constexpr double kDistAccretor = 3.88190065213158 * 1e10; ///< характерное расстояние
constexpr double kStefanBoltzmann = 5.670374419 * 1e-5;   ///< постоянная Стефана-Больцмана[ эрг·с^−1·см^−2·К^−4]
constexpr double kM_electron = 9.109383701528 * 1e-28;    ///< масса электрона [г]
constexpr double kR_electron = 2.8179403 * 1e-13;         ///< радиус электрона [см]
constexpr double kGravity = 6.6726 * 1e-8;                ///< гравитационная постоянная [м3·с−2·кг−1]

#else //СИ

constexpr double kC_Light = 299792458.0;                 ///<скорость света в м/c
constexpr double kC_LightInv = (1.0 / (kC_Light));       ///< обратная величина к скорости света
constexpr double kR_gas = 8.314;                         ///< газовая постоянная [ Дж/(моль*К)]
constexpr double kH_plank = 6.62 * 1e-34;                ///< постоянная Планка[кг * м^2 /с]
constexpr double k_boltzmann = 1.38 * 1e-23;             ///< постоянная Больцмана[Дж/K] = [ кг*м^2/(с^2*T)]
constexpr double kSigma_thomson = 6.65210 * 1e-29;       ///< сечение томсоновского рассеяния [m^2]
constexpr double kM_hydrogen = 1.6735575 * 1e-27;        ///< масса водорода[кг]
constexpr double kStefanBoltzmann = 5.670374419 * 1e-8;  ///< постоянная Стефана-Больцмана[Вт*м^-2*К^−4]
constexpr double kM_electron = 9.109383701528 * 1e-31;   ///< масса электрона [кг]
constexpr double kR_electron = 2.8179403 * 1e-15;        ///< радиус электрона [м]
constexpr double kGravity = 6.67430151515151515 * 1e-11; ///< гравитационная постоянная [м3·с−2·кг−1]
#endif

constexpr double kStefanBoltzmann4 = 4.0 * kStefanBoltzmann;

constexpr double kEarthMass = (5.9722 * 1e25); ///< масса Земли в кг
constexpr double kSunMass = (1.9891 * 1e31);   ///< масса Солнца в кг

constexpr double kDistSun = (149.6 * 10e9); ///< расстояние до Солнца в м
constexpr double kDistMoon = 400000000.;    ///<расстояние до Луны в м

#if 1
constexpr double kDist = 1e13;         ///< характерное расстояние
constexpr double kVelocity = kC_Light; ///< характерная скорость
// constexpr double kMass =(kVelocity * kVelocity * kDist) / kGravity; ///< характерная масса
constexpr double kMass = 1e6;                 ///< характерная масса
constexpr double kTime = (kDist / kVelocity); ///< характерное время

constexpr double kDensity = (kMass / (kDist * kDist * kDist));   ///< характерная плотность
constexpr double kPressure = kMass / (kDist * kTime * kTime);    ///< характерное давление
constexpr double kRadiation = (kMass / (kTime * kTime * kTime)); ///< характерное излучение

#else
constexpr double kDist = 1;
constexpr double kMass = 1;                   ///< характерная масса
constexpr double kVelocity = 1;               ///< характерная скорость
constexpr double kTime = (kDist / kVelocity); ///< характерное время
constexpr double kDensity = 1;
constexpr double kPressure = (kDensity * kVelocity * kVelocity);              ///< характерное давление
constexpr double kRadiation = (kDensity * kVelocity * kVelocity * kVelocity); ///< характерное излучение;
#endif

#endif //! GLOBAL_CONSTS_H
