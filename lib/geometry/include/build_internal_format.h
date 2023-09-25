#ifndef BUILD_INTERNAL_FORMAT_H
#define BUILD_INTERNAL_FORMAT_H

#include "prj_config.h"
#ifdef USE_VTK

#include "json_struct.h"

int BuildDataFromVTK(const global_files_t &glb_files);

#endif

#endif //! BUILD_INTERNAL_FORMAT_H