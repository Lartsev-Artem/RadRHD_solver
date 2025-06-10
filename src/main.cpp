#include "mpi_ext.h"
#include "reader_json.h"

#include "reader_bin.h"
#include "writer_bin.h"

int main(int argc, char *argv[])
{
  MPI_START(argc, argv);
  INIT_ENVIRONMENT(argc, argv);

  //...
  std::vector<Type> c;
  files_sys::bin::ReadSimple(glb_files.base_address + F_CENTERS, c);

  std::vector<Type> cnt(c.size());
  for (size_t i = 0; i < cnt.size(); i++)
  {
    cnt[i] = i;
  }
  files_sys::bin::WriteSimple(glb_files.base_address + "cnt.bin", cnt);

  DEINIT_ENVIRONMENT(argc, argv);
  MPI_END;
  return e_completion_success;
}