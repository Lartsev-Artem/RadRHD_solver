#ifndef READER_JSON
#define READER_JSON

#include <string>

#include "global_value.h"

#include "../json/json_struct.h"

namespace files_sys{
namespace json{

int ReadFileJson(const std::string& file, global_files_t& st);
int ReadFileJson(const std::string& file, solve_mode_t& st);
int ReadFileJson(const std::string& file, hllc_value_t& st);

template <typename file>
int ReadStartSettings(file file_set, global_files_t& glb_files, solve_mode_t* solve_mode = nullptr, hllc_value_t* hllc_conf = nullptr)
{		
	if (ReadFileJson(file_set, glb_files))
	{
		return e_completion_fail;
	} 
	
	glb_files.name_file_settings = file_set;
	
	glb_files.Build();

#ifdef DEBUG
	glb_files.print();
#endif

	if (solve_mode != nullptr)
	{
		if (ReadFileJson(glb_files.solve_configuration, *solve_mode)) return 1;
	}

	if (hllc_conf != nullptr)
	{
		if (ReadFileJson(glb_files.name_file_hllc_set, *hllc_conf)) return 1;
	}
	return e_completion_success;
}

} // namespace json
} //namespace files_sys
#endif // !READER_JSON
