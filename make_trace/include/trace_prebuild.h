#if !defined TRACE_PREBUILD_H && defined MAKE_TRACE
#define TRACE_PREBUILD_H

#include "json_struct.h"

namespace trace {

int PreBuild(const global_files_t &glb_files);

} // namespace trace
#endif //! TRACE_PREBUILD_H