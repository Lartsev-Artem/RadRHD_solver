#ifdef BUILD_GRAPH

#include "global_value.h"
#include "graph_config.h"

#ifdef USE_STRANGE_FUNCTION
int ReadInitBoundarySet(const std::string name_file_boundary, std::set<IntId> &boundary_cells) {

  std::ifstream ifile;

  ifile.open(name_file_boundary);
  if (!ifile.is_open()) {
    std::cout << "Error read file boundary\n";
    return 1;
  }

  IntId N;
  ifile >> N;
  // std::cout << "Boundary has " << N << "cells\n";
  while (ifile >> N) {
    boundary_cells.emplace(N);
  }

  ifile.close();
  return 0;
}
int ReadInnerBoundary(const std::string name_file_boundary_inner, std::set<IntId> &id_inner_boundary_face) {

  std::ifstream ifile;

  ifile.open(name_file_boundary_inner);
  if (!ifile.is_open()) {
    std::cout << "Error read file boundary_inner\n";
    return 1;
  }

  id_inner_boundary_face.clear();

  IntId buf;
  ifile >> buf;
  //	std::cout << "Inner boundary has " << buf << "faces\n";
  while (ifile >> buf) {
    id_inner_boundary_face.emplace(buf);
  }

  ifile.close();
  return 0;
}
int ReadInnerCellOfSphere(const std::string name_file_inner_sphere, std::vector<Face> &inner_faces) {

  std::ifstream ifile;

  ifile.open(name_file_inner_sphere);
  if (!ifile.is_open()) {
    std::cout << "Error read inner_sphere\n";
    return 1;
  }

  int N;
  ifile >> N;

  inner_faces.resize(N);

  Face face;
  for (int i = 0; i < N; ++i) {

    ifile >> face.A[0] >> face.A[1] >> face.A[2] >> face.B[0] >> face.B[1] >> face.B[2] >> face.C[0] >> face.C[1] >> face.C[2];
    inner_faces[i] = face;
  }

  ifile.close();
  return 0;
}
int ReadInnerCellOfSphereAndId(const std::string name_file_face_and_id, std::map<IntId, Face> &inner_faces) {

  std::ifstream ifile;

  ifile.open(name_file_face_and_id);
  if (!ifile.is_open()) {
    std::cout << "Error read inner_sphere\n";
    return 1;
  }

  int N;
  ifile >> N;

  Face face;
  IntId id;
  for (int i = 0; i < N; ++i) {
    ifile >> id;
    ifile >> face.A[0] >> face.A[1] >> face.A[2] >> face.B[0] >> face.B[1] >> face.B[2] >> face.C[0] >> face.C[1] >> face.C[2];
    inner_faces.emplace(id, face);
  }

  ifile.close();
  return 0;
}

int WriteFileGraph(const int i, const std::string &name_file_graph, const std::vector<IntId> &graph) {

  std::unique_ptr<FILE, int (*)(FILE *)> file_graph(fopen((name_file_graph + std::to_string(i) + ".bin").c_str(), "wb"), fclose);
  if (!file_graph) {
    printf("file_graph is not opened for writing\n");
    return 1;
  }

  const int n = graph.size();
  fwrite_unlocked(graph.data(), sizeof(IntId), n, file_graph.get());

  fclose(file_graph.get());

  std::unique_ptr<FILE, int (*)(FILE *)> file_id(fopen((std::string(glb_files.base_adress) + "id_defining_faces" + std::to_string(i) + ".bin").c_str(), "wb"), fclose);
  if (!file_id) {
    printf("file_id is not opened for writing\n");
    return 1;
  }

  int size = id_try_surface.size();
  fwrite_unlocked(&size, sizeof(int), 1, file_id.get());
  fwrite_unlocked(id_try_surface.data(), sizeof(IntId), size, file_id.get());

  fclose(file_id.get());

  std::unique_ptr<FILE, int (*)(FILE *)> file_dist(fopen((std::string(glb_files.base_adress) + "dist_defining_faces" + std::to_string(i) + ".bin").c_str(), "wb"), fclose);
  if (!file_dist) {
    printf("file_dist is not opened for writing\n");
    return 1;
  }

  size = dist_try_surface.size();
  fwrite_unlocked(&size, sizeof(int), 1, file_dist.get());
  fwrite_unlocked(dist_try_surface.data(), sizeof(Type), size, file_dist.get());

  fclose(file_dist.get());

  std::unique_ptr<FILE, int (*)(FILE *)> file_x(fopen((std::string(glb_files.base_adress) + "x_defining_faces" + std::to_string(i) + ".bin").c_str(), "wb"), fclose);
  if (!file_x) {
    printf("file_x is not opened for writing\n");
    return 1;
  }

  size = x_try_surface.size();
  fwrite_unlocked(&size, sizeof(int), 1, file_x.get());
  fwrite_unlocked(x_try_surface.data(), sizeof(Vector3), size, file_x.get());

  fclose(file_x.get());

  return 0;
}

int WriteFileGraph(std::unique_ptr<FILE, int (*)(FILE *)> &file_graph, std::unique_ptr<FILE, int (*)(FILE *)> &file_id,
                   std::unique_ptr<FILE, int (*)(FILE *)> &file_dist, std::unique_ptr<FILE, int (*)(FILE *)> &file_x,
                   const int i, const int n, const std::vector<IntId> &graph) {

  fwrite_unlocked(graph.data(), sizeof(IntId), n, file_graph.get());

  id_try_size += id_try_surface.size();
  fwrite_unlocked(id_try_surface.data(), sizeof(IntId), id_try_surface.size(), file_id.get());

  dist_try_size += dist_try_surface.size();
  fwrite_unlocked(dist_try_surface.data(), sizeof(Type), dist_try_surface.size(), file_dist.get());

  x_try_size += x_try_surface.size();
  fwrite_unlocked(x_try_surface.data(), sizeof(Vector3), x_try_surface.size(), file_x.get());

  return 0;
}
#endif

#endif //!  BUILD_GRAPH