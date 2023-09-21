/**
 * @file dbgdef.h
 * @brief Файл содержит определения глобальных отладочных макросов
 *
 */
#ifndef DBG_DEF
#define DBG_DEF

#include <stdio.h>

#include <fstream>
#include <string>

#include "global_types.h"
#include "global_value.h"
#include "prj_config.h"

#define Files_log std::string(glb_files.base_address + F_LOG).c_str()

/**
 * @brief Вывод сообщения с указанием позиции в коде
 * @param _file файл лога
 * @param _msg сообщение
 *
 */
#define WRITE_POS(_file, _msg)                                       \
  {                                                                  \
    std::ofstream out(_file, std::ios::app);                         \
    out << _msg << " :: " << __FILE__ << " " << __FUNCTION__ << ", " \
        << __LINE__ << "c.\n";                                       \
    out.close();                                                     \
  }

/**
 * @brief Печать позиции в коде
 *
 */
#define D_L WRITE_POS(Files_log, "D_L: ")

/**
 * @brief Прерывание работы программы с указанием места
 *
 */
#define D_LD WRITE_POS(Files_log, "DIE: ") abort();

/**
 * @brief Условное прерывание работы программы с указанием места
 *
 */
#define DIE_IF(_cond) \
  if (_cond)          \
    D_LD;

/**
 * @brief Лог ошибок. вывод в файл
 * \note Лог ошибок включен всегда
 */
#define WRITE_LOG_ERR(...)                       \
  {                                              \
    char buf[1024];                              \
    sprintf(buf, __VA_ARGS__);                   \
    std::ofstream out(Files_log, std::ios::app); \
    out << buf << "\n";                          \
    out.close();                                 \
  }

/**
 * @brief Лог ошибок. вывод на экран
 *
 */
#define PRINT_LOG_ERR(...)     \
  {                            \
    char buf[1024];            \
    sprintf(buf, __VA_ARGS__); \
    printf("%s\n", buf);       \
  }

/**
 * @brief Прерывание работы программы с выводом сообщения
 *
 */
#define EXIT_ERR(...)           \
  {                             \
    WRITE_LOG_ERR(__VA_ARGS__); \
    D_LD;                       \
  }

/**
 * @brief Возврат из функции с ошибкой и записью лога
 *
 */
#define RETURN_ERR(...)        \
  {                            \
    WRITE_LOG_ERR(__VA_ARGS__) \
    return e_completion_fail;  \
  }

#ifdef WRITE_GLOBAL_LOG

#ifdef WRITE_MPI_LOG
#include "mpi_ext.h"

/**
 * @brief Вывод лога в отдельный для каждого узла файл
 *
 */
#define WRITE_LOG(str)                                            \
  {                                                               \
    std::ofstream ofile;                                          \
    ofile.open(Files_log + std::to_string(get_mpi_id()) + ".txt", \
               std::ios::app);                                    \
    ofile << str;                                                 \
    ofile.close();                                                \
  }

#else

#define WRITE_LOG(str) WRITE_LOG_ERR(str)

#endif // WRITE_MPI_LOG

#else
#define CONVERT_TO_STRING(s, ...) #s #__VA_ARGS__
#define WRITE_LOG(str) \
  { CONVERT_TO_STRING(str); }
#endif // WRITE_GLOBAL_LOG

//#undef Files_log
#undef WRITE_POS
#endif // DBG_DEF