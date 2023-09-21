#ifndef GLOBAL_DEF
#define GLOBAL_DEF

#include "dbgdef.h"
#include "prj_config.h"


#define CELL_SIZE (NUMBER_OF_MEASUREMENTS + 1)

#define CONVERT_TO_STRING(s, ...) #s #__VA_ARGS__

#define OPEN_FSTREAM(file, namefile) \
  file.open(namefile);               \
  if (!file.is_open()) RETURN_ERR("Error : file %s is not open", namefile);

#define OPEN_FILE(file, namefile, mod) \
  file = fopen(namefile, mod);         \
  if (!file) RETURN_ERR("Error : file %s is not open", namefile);

#define check_bit(word, idx) (((word) >> (idx)) & 0x1)  // проверка i-го бита
#define clear_bit(word, idx) ((word) & (~(1 << (idx))))  // выключение i-го бита
#define set_bit(word, idx) ((word) | (1 << (idx)))  // установка i-го бита

#define SIGN(a) (a < 0.0 ? -1.0 : 1.0)

#define CREATE_STRUCT(name, type, ...) \
  struct name {                        \
    type __VA_ARGS__;                  \
  }

#define PRINT_STRUCT(st, type)                           \
  type* str = (type*)&st;                                \
  while (str < (type*)&st + sizeof(st) / sizeof(type)) { \
    std::cout << *str++ << '\n';                         \
  }

#define FILL_STRUCT(st, type, val)                         \
  {                                                        \
    type* str = (type*)&st;                                \
    while (str < (type*)&st + sizeof(st) / sizeof(type)) { \
      *str++ = val;                                        \
    }                                                      \
  }

///\todo: file module

#define WRITE_FILE_VECTOR(name_file, data, value)         \
  {                                                       \
    FILE* f;                                              \
    int n = data.size();                                  \
    f = fopen(name_file, "wb");                           \
    if (!f) RETURN_ERRS("file %s not open\n", name_file); \
    fwrite(&n, sizeof(int), 1, f);                        \
    for (auto& el : data) {                               \
      fwrite(&el.value, sizeof(el.value), 1, f);          \
    }                                                     \
    fclose(f);                                            \
  }

#define WRITE_FILE(name_file, data, n)                    \
  {                                                       \
    FILE* f;                                              \
    f = fopen(name_file, "wb");                           \
    if (!f) RETURN_ERRS("file %s not open\n", name_file); \
    fwrite(&n, sizeof(int), 1, f);                        \
    fwrite(data, sizeof(data[0]), n, f);                  \
    fclose(f);                                            \
  }

#define WRITE_FILE_PHYS(name_file, data, value, type, param) \
  {                                                          \
    FILE* f;                                                 \
    int n = data.size();                                     \
    f = fopen(name_file, "wb");                              \
    if (!f) RETURN_ERRS("file %s not open\n", name_file);    \
    fwrite(&n, sizeof(int), 1, f);                           \
    for (auto& el : data) {                                  \
      type x = el.value * param;                             \
      fwrite(&x, sizeof(x), 1, f);                           \
    }                                                        \
    fclose(f);                                               \
  }

#endif