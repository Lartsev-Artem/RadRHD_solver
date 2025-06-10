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
#include "prj_config.h"

#include "json_struct.h"
#include "mpi_ext.h"
#include "timer.h"
extern global_files_t glb_files;
#define Files_log std::string(glb_files.base_address + "File_Logs.txt").c_str()

#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

/**
 * @brief Вывод сообщения с указанием позиции в коде
 * @param _file файл лога
 * @param _msg сообщение
 *
 */
#define WRITE_POS(_file, _msg)                                                            \
  {                                                                                       \
    std::ofstream out(_file, std::ios::app);                                              \
    out << Timer::get_time() << _msg << " :: " << __FILE__ << " " << __FUNCTION__ << ", " \
        << __LINE__ << "c.\n";                                                            \
    out.close();                                                                          \
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
#define D_LD                      \
  {                               \
    WRITE_POS(Files_log, "DIE: ") \
    abort();                      \
  }

/**
 * @brief Условное прерывание работы программы с указанием места
 *
 */
#define DIE_IF(_cond)             \
  if (UNLIKELY(_cond))            \
  {                               \
    WRITE_POS(Files_log, #_cond); \
    D_LD;                         \
  }

#ifdef DEBUG
#define STOP_IF(_cond)            \
  if (UNLIKELY(_cond))            \
  {                               \
    WRITE_POS(Files_log, #_cond); \
    D_LD;                         \
  }
#else
#define STOP_IF(_cond)
#endif

#define DIE_IF_ACTION(_cond, _act) \
  if (UNLIKELY(_cond))             \
  {                                \
    _act;                          \
    D_LD;                          \
  }

#ifdef LOG_OUT_TO_SCREEN
#define WRITE_LOG_ERR(...) printf(__VA_ARGS__);
#else
/**
 * @brief Лог ошибок. вывод в файл
 * \note Лог ошибок включен всегда
 */
#define WRITE_LOG_ERR(...)                                                                     \
  {                                                                                            \
    char buf[1024];                                                                            \
    sprintf(buf, __VA_ARGS__);                                                                 \
    std::ofstream out(Files_log, std::ios::app);                                               \
    out << Timer::get_time() << "from " << std::to_string(get_mpi_id()) + ": " << buf << "\n"; \
    out.close();                                                                               \
  }
#endif

/**
 * @brief Лог ошибок. вывод на экран
 *
 */
#define PRINT_LOG_ERR(...)                             \
  {                                                    \
    char buf[1024];                                    \
    sprintf(buf, __VA_ARGS__);                         \
    printf("%s %s\n", Timer::get_time().c_str(), buf); \
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
// #include "mpi_ext.h"

/**
 * @brief Вывод лога в отдельный для каждого узла файл
 *
 */
#define WRITE_LOG(...)                                            \
  {                                                               \
    char buf[1024];                                               \
    sprintf(buf, __VA_ARGS__);                                    \
    std::ofstream ofile;                                          \
    ofile.open(Files_log + std::to_string(get_mpi_id()) + ".txt", \
               std::ios::app);                                    \
    ofile << Timer::get_time() << buf;                            \
    ofile.close();                                                \
  }

#else

#define WRITE_LOG(...) WRITE_LOG_ERR(__VA_ARGS__)

#endif // WRITE_MPI_LOG

#else
#define CONVERT_TO_STRING(s, ...) #s #__VA_ARGS__
#define WRITE_LOG(...)              \
  {                                 \
    CONVERT_TO_STRING(__VA_ARGS__); \
  }
#endif // WRITE_GLOBAL_LOG

#ifdef DEBUG
#include <unistd.h>
#define GDB_ATTACH      \
  {                     \
    volatile int i = 0; \
    while (i == 0)      \
      sleep(1);         \
  }

#define WRITE_CELL_INFO(id, _grid)                                                                                                                 \
  {                                                                                                                                                \
    geo_cell_t *cell = &_grid.cells[id].geo;                                                                                                       \
    geo_face_t *f1 = &_grid.faces[cell->id_faces[0]].geo;                                                                                          \
    geo_face_t *f2 = &_grid.faces[cell->id_faces[1]].geo;                                                                                          \
    geo_face_t *f3 = &_grid.faces[cell->id_faces[2]].geo;                                                                                          \
    geo_face_t *f4 = &_grid.faces[cell->id_faces[3]].geo;                                                                                          \
    WRITE_LOG("N_bits=%u, V=%lf, id_faces: %d %d %d %d\n", cell->sign_n.bits, cell->V, cell->id_faces[0],                                          \
              cell->id_faces[1], cell->id_faces[2], cell->id_faces[3]);                                                                            \
    WRITE_LOG("f1: n=[%lf, %lf, %lf], S=%lf, neigh: %d, %d is_reg=%d\n", f1->n[0], f1->n[1], f1->n[2], f1->S, f1->id_l, f1->id_r, f1->is_regular); \
    WRITE_LOG("f2: n=[%lf, %lf, %lf], S=%lf, neigh: %d, %d is_reg=%d\n", f2->n[0], f2->n[1], f2->n[2], f2->S, f2->id_l, f2->id_r, f2->is_regular); \
    WRITE_LOG("f3: n=[%lf, %lf, %lf], S=%lf, neigh: %d, %d is_reg=%d\n", f3->n[0], f3->n[1], f3->n[2], f3->S, f3->id_l, f3->id_r, f3->is_regular); \
    WRITE_LOG("f4: n=[%lf, %lf, %lf], S=%lf, neigh: %d, %d is_reg=%d\n", f4->n[0], f4->n[1], f4->n[2], f4->S, f4->id_l, f4->id_r, f4->is_regular); \
  }

#endif

// #undef Files_log
// #undef WRITE_POS
#endif // DBG_DEF