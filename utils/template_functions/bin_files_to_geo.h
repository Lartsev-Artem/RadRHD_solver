#ifndef BIN_FILES_TO_GEO
#define BIN_FILES_TO_GEO

#include "convert_bin_to_geo.h"
#include "global_value.h"
#include "reader_bin.h"
#include "writer_bin.h"

template <typename file>
int BinToGeo(file base_address) {

  std::vector<IntId> neighbours_id_faces;
  std::vector<Normals> normals;
  std::vector<Type> areas_faces;
  std::vector<Type> volume;
  std::vector<Vector3> centers;

  if (files_sys::bin::ReadSimple(base_address + F_NEIGHBOR, neighbours_id_faces))
    RETURN_ERR("Error reading file neighbours\n");
  if (files_sys::bin::ReadNormals(base_address + F_NORMALS, normals))
    RETURN_ERR("Error reading file normals\n");
  if (files_sys::bin::ReadSimple(base_address + F_AREAS, areas_faces))
    RETURN_ERR("Error reading file areas_faces\n");
  if (files_sys::bin::ReadSimple(base_address + F_VOLUME, volume))
    RETURN_ERR("Error reading file volume\n");
  if (files_sys::bin::ReadSimple(base_address + F_CENTERS, centers))
    RETURN_ERR("Error reading file centers\n");

  std::vector<elem_t> cells;
  std::vector<face_t> faces;
  ConvertBinToGeo(neighbours_id_faces, normals, areas_faces, volume, centers, faces, cells);

  uint32_t err = files_sys::bin::WriteGridGeo(base_address + F_GEO_CELLS, cells);
  err |= files_sys::bin::WriteGridGeo(base_address + F_GEO_FACES, faces);
  return err;
}
#endif //!