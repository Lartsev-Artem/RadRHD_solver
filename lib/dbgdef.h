#ifndef DBG_DEF
#define DBG_DEF

#include <stdio.h>

#include <fstream>
#include <string>

#include "global_types.h"
#include "prj_config.h"

// todo: требует glob_value.h!!!
#include "global_value.h"
#define Files_log std::string(glb_files.base_address + F_LOG).c_str()

#define WRITE_POS(_file, _msg)                                       \
  {                                                                  \
    std::ofstream out(_file, std::ios::app);                         \
    out << _msg << " :: " << __FILE__ << " " << __FUNCTION__ << ", " \
        << __LINE__ << "c.\n";                                       \
    out.close();                                                     \
  }

#define D_L WRITE_POS(Files_log, "D_L: ")
#define D_LD WRITE_POS(Files_log, "DIE: ") abort();
#define DIE_IF(_cond) \
  if (_cond) D_LD;

// Лог ошибок включен всегда
#define WRITE_LOG_ERR(...)                       \
  {                                              \
    char buf[1024];                              \
    sprintf(buf, __VA_ARGS__);                   \
    std::ofstream out(Files_log, std::ios::app); \
    \  
    out << buf                                   \
        << "\n";                                 \
    out.close();                                 \
  }

#define PRINT_LOG_ERR(...)     \
  {                            \
    char buf[1024];            \
    sprintf(buf, __VA_ARGS__); \
    printf("%s\n", buf);       \
  }

#define EXIT_ERR(...)           \
  {                             \
    WRITE_LOG_ERR(__VA_ARGS__); \
    D_LD;                       \
  }
#define RETURN_ERR(...) \
  { WRITE_LOG_ERR(__VA_ARGS__) return e_completion_fail; }

#ifdef WRITE_GLOBAL_LOG

#ifdef WRITE_MPI_LOG
#include "mpi_ext.h"

#define WRITE_LOG(str)                                            \
  {                                                               \
    std::ofstream ofile;                                          \
    ofile.open(Files_log + std::to_string(get_mpi_id()) + ".txt", \
               std::ios::app);                                    \
    ofile << str;                                                 \
    ofile.close();                                                \
  }

#define DBG_POINT                                                            \
  WRITE_LOG("dbg_p: " << __FILE__ << " " << __FUNCTION__ << ", " << __LINE__ \
                      << "c.\n")

#else

#define WRITE_LOG(str) WRITE_LOG_ERR(str)
#define DBG_POINT                                                             \
  WRITE_LOG("dbg_p: ", << __FILE__ << " " << __FUNCTION__ << ", " << __LINE__ \
                       << "c.\n")

#endif  // WRITE_MPI_LOG

#else
#define CONVERT_TO_STRING(s, ...) #s #__VA_ARGS__
#define WRITE_LOG(str) \
  { CONVERT_TO_STRING(str); }
#endif  // WRITE_GLOBAL_LOG

//#undef Files_log
#undef WRITE_POS
#endif  // DBG_DEF