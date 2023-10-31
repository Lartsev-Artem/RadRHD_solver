#include "reader_txt.h"
#include "mpi_ext.h"
#include "mpi_shifts.h"
#include <map>

int files_sys::txt::ReadSphereDirectionCartesian(const std::string &file_sphere_direction,
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

  std::vector<IdType> send;
  std::vector<IdType> disp;

  GetDisp(get_mpi_np(), N, disp);
  GetSend(get_mpi_np(), N, send);

  grid_direction.loc_size = send[get_mpi_id()];
  grid_direction.loc_shift = disp[get_mpi_id()];

  // WRITE_LOG("Read %s success\n", file_sphere_direction.c_str());
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

int files_sys::txt::ReadTableFunc(const std::string &file_func, TableFunc &tab_func) {

  std::ifstream ifile;
  OPEN_FSTREAM(ifile, file_func.c_str());

  std::string str;

  /*--------------Read preamble-------------*/

  std::getline(ifile, str);
  int Nd = std::stoi(&str[str.find("=", 0)] + 1);

  std::getline(ifile, str);
  int Nt = std::stoi(&str[str.find("=", 0)] + 1);

  std::getline(ifile, str);
  std::getline(ifile, str);

  TableFunc lambda(Nd, Nt);

  std::vector<Type> density(Nd);
  std::vector<Type> temperature(Nt);

  /*--------------Read data-------------*/

  for (int i = 0; i < Nd; i++) {
    ifile >> density[i] >> temperature[0] >> lambda.data[i * Nt + 0];
    for (int j = 1; j < Nt; j++) {
      Type buf;
      ifile >> buf >> temperature[j] >> lambda.data[i * Nt + j];
    }
  }
  ifile.close();

  lambda.step_x = density[1] - density[0];
  lambda.step_y = temperature[1] - temperature[0];

  lambda.max_x = -std::numeric_limits<Type>::max();
  lambda.min_x = std::numeric_limits<Type>::max();
  for (int i = 0; i < Nd; i++) {
    lambda.max_x = std::max(lambda.max_x, density[i]);
    lambda.min_x = std::min(lambda.min_x, density[i]);
  }

  lambda.max_y = -std::numeric_limits<Type>::max();
  lambda.min_y = std::numeric_limits<Type>::max();

  for (int i = 0; i < Nt; i++) {
    lambda.max_y = std::max(lambda.max_y, temperature[i]);
    lambda.min_y = std::min(lambda.min_y, temperature[i]);
  }

  WRITE_LOG("Tab func: %dx%d: hx=%lf hy=%lf maxX= %lf maxY=%lf minX= %lf minY= %lf\n",
            lambda.Nx, lambda.Ny,
            lambda.step_x, lambda.step_y,
            lambda.max_x, lambda.max_y,
            lambda.min_x, lambda.min_y);

  return e_completion_success;
}