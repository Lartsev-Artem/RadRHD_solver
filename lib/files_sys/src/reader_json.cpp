#include <fstream>
#include <string>

#include "global_def.h"
#include "../json/json_struct.h"

namespace files_sys{
namespace json{

static int ReadFileJson(nlohmann::json& j, const std::string& file)
{
	std::ifstream ifile(file);
	OPEN_FSTREAM(ifile, file.c_str());

	ifile >> j;
	ifile.close();
	return 0;
}

#define READ_JSON \
nlohmann::json j; \
ReadFileJson(j, file); \
if (j.is_null()) \
{ \
	RETURN_ERR("Error : j.is_null\n"); \
} \
st = j; \
return e_completion_success; \

int Read(const std::string& file, global_files_t& st)
{
	READ_JSON;
}
int Read(const std::string& file, solve_mode_t& st)
{
	READ_JSON;
}
int Read(const std::string& file, hllc_value_t& st)
{
	READ_JSON;
}
} // namespace json
} //namespace files_sys