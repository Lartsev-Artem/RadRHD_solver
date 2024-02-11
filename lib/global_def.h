/**
 * @file global_def.h
 * @brief Файл содержит глобальные макросы
 *
 */
#ifndef GLOBAL_DEF
#define GLOBAL_DEF

#include "dbgdef.h"
#include "prj_config.h"

/**
 * @brief Число граней элемента
 *
 */
#define CELL_SIZE (NUMBER_OF_MEASUREMENTS + 1)

#define NODE_SIZE 3 ///< число узлов интерполяции на грани

#define CONVERT_TO_STRING(s, ...) #s #__VA_ARGS__

/**
 * @brief Безопасное открытие файлового потока
 *
 */
#define OPEN_FSTREAM(file, namefile) \
  file.open(namefile);               \
  if (UNLIKELY(!file.is_open()))     \
    RETURN_ERR("Error : file %s is not open", namefile);

/**
 * @brief Безопасное открытие файла
 *
 */
#define OPEN_FILE(file, namefile, mod) \
  file = fopen(namefile, mod);         \
  if (UNLIKELY(!file))                 \
    RETURN_ERR("Error : file %s is not open", namefile);

/**
 * @brief проверка i-го бита
 *
 */
#define CHECK_BIT(word, idx) (((word) >> (idx)) & 0x1)

/**
 * @brief выключение i-го бита
 *
 */
#define CLEAR_BIT(word, idx) ((word) & (~(1 << (idx))))

/**
 * @brief  установка i-го бита
 *
 */
#define SET_BIT(word, idx) ((word) | (1 << (idx)))

/**
 * @brief Знак величины
 *
 */
#define SIGN(a) (a < 0.0 ? -1.0 : 1.0)

/**
 * @brief Быстрое создание структуры с элементами одного типа
 *
 */
#define CREATE_STRUCT(name, type, ...) \
  struct name {                        \
    type __VA_ARGS__;                  \
  }

/**
 * @brief Быстрая печать структуры с элементами одного типа
 *
 */
#define PRINT_STRUCT(st, type)                            \
  type *str = (type *)&st;                                \
  while (str < (type *)&st + sizeof(st) / sizeof(type)) { \
    std::cout << *str++ << '\n';                          \
  }

/**
 * @brief Быстрое заполнение структуры с элементами одного типа
 *
 */
#define FILL_STRUCT(st, type, val)                          \
  {                                                         \
    type *str = (type *)&st;                                \
    while (str < (type *)&st + sizeof(st) / sizeof(type)) { \
      *str++ = val;                                         \
    }                                                       \
  }

#endif