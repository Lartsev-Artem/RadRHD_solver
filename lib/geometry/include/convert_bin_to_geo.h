/**
 * @file convert_bin_to_geo.h
 * @brief Функция перевода бинарных файлов геометрии в структуру сетки
 *
 */
#ifndef CONVERT_BIN_TO_GEO_H
#define CONVERT_BIN_TO_GEO_H

#include "geo_types.h"
#include "solvers_struct.h"

void ConvertBinToGeo(std::vector<IntId> &neighbours_id_faces,
                     std::vector<Normals> &normals, std::vector<Type> &areas_faces,
                     std::vector<Type> &volume, std::vector<Vector3> &centers,
                     std::vector<face_t> &faces, std::vector<elem_t> &cells);

#endif //! CONVERT_BIN_TO_GEO_H