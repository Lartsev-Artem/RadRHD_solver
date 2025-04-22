#include "reader_json.h"

#include <fstream>
#include <string>

#include "global_def.h"
#include "json_struct.h"

static int ReadFileJson(nlohmann::json &j, const std::string &file) {
  std::ifstream ifile(file);
  OPEN_FSTREAM(ifile, file.c_str());

  ifile >> j;
  ifile.close();
  return 0;
}

#define READ_JSON                      \
  nlohmann::json j;                    \
  ReadFileJson(j, file);               \
  if (j.is_null()) {                   \
    RETURN_ERR("Error : j.is_null\n"); \
  }                                    \
  st = j;                              \
  return e_completion_success;

int files_sys::json::Read(const std::string &file, global_files_t &st) { READ_JSON; }
int files_sys::json::Read(const std::string &file, solve_mode_t &st) { READ_JSON; }
int files_sys::json::Read(const std::string &file, hllc_value_t &st) { READ_JSON; }
