#include "reader_txt.h"
#include <map>

int files_sys::txt::ReadSphereDirectionСartesian(const std::string &file_sphere_direction,
                                                 grid_directions_t &grid_direction) {
  std::ifstream ifile;
  OPEN_FSTREAM(ifile, file_sphere_direction.c_str());

  int N = 0;
  ifile >> N;

  grid_direction.size = N;
  grid_direction.directions.resize(N);

  for (int i = 0; i < N; i++) {
    ifile >> grid_direction.directions[i].area;
    ifile >> grid_direction.directions[i].dir[0] >>
        grid_direction.directions[i].dir[1] >>
        grid_direction.directions[i].dir[2];
  }
  ifile >> grid_direction.full_area;
  ifile.close();

  return e_completion_success;
}

int files_sys::txt::ReadInitBoundarySetInFaces(const std::string &file_face_id, std::map<IntId, FaceCell> &inter_faces) {

  std::ifstream ifile;
  OPEN_FSTREAM(ifile, file_face_id.c_str());

  int N;
  ifile >> N;

  inter_faces.clear();

  FaceCell face;
  for (int i = 0; i < N; ++i) {
    ifile >> face;
    inter_faces.emplace(face.face_id / 4, face);
  }

  ifile.close();
  return e_completion_success;
}