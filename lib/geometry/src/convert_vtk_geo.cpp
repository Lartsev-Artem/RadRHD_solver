#include "convert_vtk_geo.h"
#ifdef USE_VTK
#include "global_value.h"

#include <vtkCell.h>
#include <vtkIdList.h>
#include <vtkPoints.h>

int GetFacesPoints(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Face> &faces) {

  const vtkIdType n = unstructured_grid->GetNumberOfCells();
  faces.resize(n * CELL_SIZE);

  for (vtkIdType i = 0; i < n; i++) {
    for (vtkIdType j = 0; j < CELL_SIZE; j++) {
      for (vtkIdType k = 0; k < CELL_SIZE - 1; k++)
        faces[i * CELL_SIZE + j][k] = Vector3(unstructured_grid->GetCell(i)->GetFace(j)->GetPoints()->GetPoint(k));
    }
  }

  return e_completion_success;
}

int GetBoundaryCells(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::set<IntId> &boundary_cells) {

  boundary_cells.clear();
  vtkIdType N = unstructured_grid->GetNumberOfCells();

  vtkSmartPointer<vtkIdList> idc = vtkSmartPointer<vtkIdList>::New();

  for (vtkIdType i = 0; i < N; ++i) {

    for (vtkIdType j = 0; j < CELL_SIZE; ++j) {
      unstructured_grid->GetCellNeighbors(i, unstructured_grid->GetCell(i)->GetFace(j)->GetPointIds(), idc);
      if (idc->GetNumberOfIds() == 0) {
        boundary_cells.emplace((IntId)i);
        break;
      } else if (idc->GetNumberOfIds() > 1)
        RETURN_ERR("More than 1 neighbor????\n");
    }
  }

  return e_completion_success;
}

int GetInterBoundaryFacesOfSphere(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::set<IntId> &inter_boundary_faces) {

  inter_boundary_faces.clear();
  int N = unstructured_grid->GetNumberOfCells();

  vtkSmartPointer<vtkIdList> idc = vtkSmartPointer<vtkIdList>::New();

  for (int i = 0; i < N; ++i) {

    for (int j = 0; j < CELL_SIZE; ++j) {
      unstructured_grid->GetCellNeighbors(i, unstructured_grid->GetCell(i)->GetFace(j)->GetPointIds(), idc);

      if (idc->GetNumberOfIds() == 0) {

        Vector3 P(unstructured_grid->GetCell(i)->GetFace(j)->GetPoints()->GetPoint(0));

        // внешняя сфера
        if ((P - kCenterPoint).norm() < kInternalRadius) {
          inter_boundary_faces.emplace(i * CELL_SIZE + j);
        }
        break;
      } else if (idc->GetNumberOfIds() > 1)
        RETURN_ERR("More than 1 neighbor????\n");
    }
  }

  return e_completion_success;
}

int GetBoundaryFacesId(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<IntId> &boundary_faces) {

  int N = unstructured_grid->GetNumberOfCells();

  boundary_faces.clear();
  boundary_faces.reserve(N);

  vtkSmartPointer<vtkIdList> idc = vtkSmartPointer<vtkIdList>::New();

  for (int i = 0; i < N; ++i) {

    for (int j = 0; j < CELL_SIZE; ++j) {
      unstructured_grid->GetCellNeighbors(i, unstructured_grid->GetCell(i)->GetFace(j)->GetPointIds(), idc);

      if (idc->GetNumberOfIds() == 0) {
        boundary_faces.push_back(i * CELL_SIZE + j);
        break;
      } else if (idc->GetNumberOfIds() > 1)
        RETURN_ERR("More than 1 neighbor????\n");
    }
  }

  boundary_faces.shrink_to_fit();
  return e_completion_success;
}

#if NUMBER_OF_MEASUREMENTS == 3

int GetNeighborFace3D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<int> &neighbors) {

  auto GetNumberNeighborFace{[](const int a, const int b, const int c, vtkCell *neighbor_cell) {
    vtkIdList *idc;

    vtkIdType x, y, z;
    for (int i = 0; i < CELL_SIZE; i++) {
      idc = neighbor_cell->GetFace(i)->GetPointIds();
      x = idc->GetId(0);
      y = idc->GetId(1);
      z = idc->GetId(2);

      if (a == x && b == y && c == z)
        return i;
      else if (a == x && b == z && c == y)
        return i;
      else if (a == y && b == x && c == z)
        return i;
      else if (a == y && b == z && c == x)
        return i;
      else if (a == z && b == x && c == y)
        return i;
      else if (a == z && b == y && c == x)
        return i;
    }
    return (int)e_neigh_code_undef;
  }};

  int count_unique_face = 0;

  const int N = unstructured_grid->GetNumberOfCells();
  neighbors.resize(N * CELL_SIZE, e_neigh_code_undef);

  vtkSmartPointer<vtkIdList> idp = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> idc = vtkSmartPointer<vtkIdList>::New();

  int id_a, id_b, id_c;
  for (vtkIdType num_cell = 0; num_cell < N; ++num_cell) {

    for (int num_face = 0; num_face < CELL_SIZE; ++num_face) {
      if (neighbors[num_cell * CELL_SIZE + num_face] != e_neigh_code_undef)
        continue;
      ++count_unique_face;

      idp = unstructured_grid->GetCell(num_cell)->GetFace(num_face)->GetPointIds();
      id_a = idp->GetId(0);
      id_b = idp->GetId(1);
      id_c = idp->GetId(2);

      /*Может быть проблема с указателями на списки!*/
      unstructured_grid->GetCellNeighbors(num_cell, idp, idc);
      int face = num_cell * CELL_SIZE + num_face;

      if (idc->GetNumberOfIds() == 1) {
        int id_neighbor_cell = idc->GetId(0);
        int id_neighbor_face = GetNumberNeighborFace(id_a, id_b, id_c, unstructured_grid->GetCell(id_neighbor_cell));

        if (id_neighbor_face == e_neigh_code_undef) {
          RETURN_ERR("neighbor %d not found\n", (int)num_cell);
        }

        neighbors[face] = id_neighbor_cell * CELL_SIZE + id_neighbor_face;
        neighbors[id_neighbor_cell * CELL_SIZE + id_neighbor_face] = face;
      } else if (idc->GetNumberOfIds() == 0)
        neighbors[face] = e_neigh_code_out_bound; // граничная ячейка
      else {
        RETURN_ERR("More than 1 neighbor????\n");
      }
    }
  }

  return e_completion_success;
}

static void NormalAndSquareFace3D(const int cell_number, const int face_number, const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, Type &S, Vector3 &n) {

  vtkSmartPointer<vtkIdList> idp = unstructured_grid->GetCell(cell_number)->GetFace(face_number)->GetPointIds();

  double P0[3], P1[3], P2[3], P3[3]; // узлы ячейки

  unstructured_grid->GetPoint(idp->GetId(0), P0);
  unstructured_grid->GetPoint(idp->GetId(1), P1);
  unstructured_grid->GetPoint(idp->GetId(2), P2);

  double a[3], b[3]; // образующие вектора
  for (int i = 0; i < 3; i++) {
    a[i] = P1[i] - P0[i];
    b[i] = P2[i] - P0[i];
  }

  // векторное произведение -> нормаль к плоскости
  n[0] = a[1] * b[2] - a[2] * b[1];
  n[1] = -a[0] * b[2] + a[2] * b[0];
  n[2] = a[0] * b[1] - a[1] * b[0];

  S = 0.5 * n.norm(); // св-во векторного произведения
  n.normalize();      // нормировка

  vtkSmartPointer<vtkIdList> idp2 = unstructured_grid->GetCell(cell_number)->GetPointIds();

  int id; // номер узла не ячейки не принадлежащей рассматриваемой гране
  for (int i = 0; i < CELL_SIZE; i++) {
    int count = 0;
    for (int j = 0; j < 3; j++)
      if (idp2->GetId(i) != idp->GetId(j))
        count++;
    if (count == 3) {
      id = i;
      break;
    }
  }

  //Далее определение ориентации нормали

  double sum = 0;
  unstructured_grid->GetPoint(idp2->GetId(id), P3);
  /*for (size_t i = 0; i < 3; i++){
          sum += n[i] * (P3[i] - P0[i]);
  }*/

  sum = P1[0] * (P2[1] - P3[1]) * P0[2] + P0[0] * (P3[1] - P2[1]) * P1[2] +
        P0[0] * (P1[1] - P3[1]) * P2[2] + P2[2] * (P1[0] * P3[1] - P1[0] * P0[1]) +
        P3[0] * (P0[2] * (P1[1] - P2[1]) + P1[2] * (P2[1] - P0[1]) + P2[2] * (P0[1] - P1[1])) + P3[2] * (P1[0] * (P0[1] - P2[1]) + P0[0] * (P2[1] - P1[1])) +
        P2[0] * (P0[2] * (P3[1] - P1[1]) + P1[2] * (P0[1] - P3[1]) + P3[2] * (P1[1] - P0[1]));

  // разворот нормали
  if (sum < 0)
    for (int i = 0; i < 3; i++)
      n[i] = -n[i];
}

int GetNormalAndAreas3D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Normals> &normals, std::vector<Type> &areas) {

  const int n = unstructured_grid->GetNumberOfCells();
  normals.resize(n);
  areas.resize(n * (CELL_SIZE));

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < CELL_SIZE; j++) {
      NormalAndSquareFace3D(i, j, unstructured_grid, areas[i * CELL_SIZE + j], normals[i].n[j]);
    }
  }

  return e_completion_success;
}

int GetVolume3D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Type> &volumes) {

  auto GetCellVolume{[](int cell_number, const vtkSmartPointer<vtkUnstructuredGrid> &u_grid) {
    double P0[3], P1[3], P2[3], P3[3];
    vtkSmartPointer<vtkIdList> idp = u_grid->GetCell(cell_number)->GetPointIds();
    u_grid->GetPoint(idp->GetId(0), P0);
    u_grid->GetPoint(idp->GetId(1), P1);
    u_grid->GetPoint(idp->GetId(2), P2);
    u_grid->GetPoint(idp->GetId(3), P3);

    Eigen::Vector3d a, b, c; // вектора образующие тетраэдр
    for (int i = 0; i < 3; i++) {
      a[i] = P1[i] - P0[i];
      b[i] = P2[i] - P0[i];
      c[i] = P3[i] - P0[i];
    }
    return -(a.dot(b.cross(c))) / 6.0; //  св-во смешанного произведения

    /*Type a[3], b[3], c[3];
    for (size_t i = 0; i < 3; i++) {
            a[i] = P1[i] - P0[i];
            b[i] = P2[i] - P0[i];
            c[i] = P3[i] - P0[i];
    }

    Type V = a[0] * (b[1] * c[2] - c[1] * b[2]) - a[1] * (b[0] * c[2] - c[0] * b[2]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
    return fabs(V) / 6;*/
  }};

  const int n = unstructured_grid->GetNumberOfCells();
  volumes.resize(n);

  for (int i = 0; i < n; i++) {
    volumes[i] = GetCellVolume(i, unstructured_grid);
  }

  return e_completion_success;
}

static void CenterOfTetra3D(const int cell_number, const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, Vector3 &point_in_tetra) {

  // площадь грани
  auto GetArea{[](double *P0, double *P1, double *P2) {
    Eigen::Vector3d a, b;
    for (int i = 0; i < 3; i++) {
      a[i] = P1[i] - P0[i];
      b[i] = P2[i] - P0[i];
    }
    return 0.5 * a.cross(b).norm();

    // double Sum = pow(a[1] * b[2] - a[2] * b[1], 2) + pow(a[0] * b[2] - a[2] * b[0], 2) + pow(a[0] * b[1] - a[1] * b[0], 2);
    // return 0.5 * sqrt(Sum);
  }};

  double P0[3], P1[3], P2[3], P3[3];
  vtkSmartPointer<vtkIdList> idp = unstructured_grid->GetCell(cell_number)->GetPointIds();

  unstructured_grid->GetPoint(idp->GetId(0), P0);
  unstructured_grid->GetPoint(idp->GetId(1), P1);
  unstructured_grid->GetPoint(idp->GetId(2), P2);
  unstructured_grid->GetPoint(idp->GetId(3), P3);

  double Squr[4] = {GetArea(P1, P2, P3), GetArea(P0, P2, P3), GetArea(P0, P1, P3), GetArea(P0, P1, P2)};

  double Sum = Squr[0] + Squr[1] + Squr[2] + Squr[3];
  for (int i = 0; i < 3; i++) {
    point_in_tetra[i] = (Squr[0] * P0[i] + Squr[1] * P1[i] + Squr[2] * P2[i] + Squr[3] * P3[i]) / Sum;
  }
}

int GetCentersOfCells3D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Vector3> &centers) {
  const int n = unstructured_grid->GetNumberOfCells();
  centers.resize(n);

  for (int i = 0; i < n; i++) {
    CenterOfTetra3D(i, unstructured_grid, centers[i]);
  }
  return e_completion_success;
}

int GetCentersOfFaces3D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Vector3> &centers) {
  int n = unstructured_grid->GetNumberOfCells();
  centers.resize(n * CELL_SIZE);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < CELL_SIZE; j++) {
      vtkSmartPointer<vtkIdList> idp = unstructured_grid->GetCell(i)->GetFace(j)->GetPointIds();

      double P0[3], P1[3], P2[3];
      unstructured_grid->GetPoint(idp->GetId(0), P0);
      unstructured_grid->GetPoint(idp->GetId(1), P1);
      unstructured_grid->GetPoint(idp->GetId(2), P2);

      for (int k = 0; k < 3; k++) {
        centers[i * CELL_SIZE + j][k] = (P0[k] + P1[k] + P2[k]) / 3;
      }
    }
  }
  return e_completion_success;
}

int GetVertexMatrix(const size_t number_cell, const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, Eigen::Matrix4d &vertex_tetra) {

  vtkPoints *points = unstructured_grid->GetCell(number_cell)->GetPoints();
  for (int j = 0; j < 3; j++) {
    vertex_tetra(j, 0) = points->GetPoint(0)[j];
    vertex_tetra(j, 1) = points->GetPoint(1)[j];
    vertex_tetra(j, 2) = points->GetPoint(2)[j];
    vertex_tetra(j, 3) = points->GetPoint(3)[j];
  }

  for (int i = 0; i < 4; i++)
    vertex_tetra(3, i) = 1;

  return e_completion_success;
}

#elif NUMBER_OF_MEASUREMENTS == 2

int GetNeighborFace2D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<int> &neighbors) {

  auto GetNumberNeighborFace{[](const int a, const int b, vtkCell *neighbor_cell) {
    vtkIdList *idc;

    int x, y;
    for (int i = 0; i < CELL_SIZE; i++) {
      idc = neighbor_cell->GetEdge(i)->GetPointIds();
      x = idc->GetId(0);
      y = idc->GetId(1);

      if (a == x && b == y)
        return i;
      else if (a == y && b == x)
        return i;
    }
    return (int)e_neigh_code_undef;
  }};

  int count_unique_face = 0;
  const int N = unstructured_grid->GetNumberOfCells();
  neighbors.resize(N * CELL_SIZE, e_neigh_code_undef);

  vtkSmartPointer<vtkIdList> idp = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> idc = vtkSmartPointer<vtkIdList>::New();

  int id_a, id_b;
  for (vtkIdType num_cell = 0; num_cell < N; ++num_cell) {

    for (int num_face = 0; num_face < CELL_SIZE; ++num_face) {
      if (neighbors[num_cell * CELL_SIZE + num_face] != e_neigh_code_undef)
        continue;
      ++count_unique_face;

      idp = unstructured_grid->GetCell(num_cell)->GetEdge(num_face)->GetPointIds();
      id_a = idp->GetId(0);
      id_b = idp->GetId(1);

      /*Может быть проблема с указателями на списки!*/
      unstructured_grid->GetCellNeighbors(num_cell, idp, idc);
      int face = num_cell * CELL_SIZE + num_face;

      if (idc->GetNumberOfIds() == 1) {
        int id_neighbor_cell = idc->GetId(0);
        int id_neighbor_face = GetNumberNeighborFace(id_a, id_b, unstructured_grid->GetCell(id_neighbor_cell));

        if (id_neighbor_face == e_neigh_code_undef) {
          RETURN_ERR("neighbor %d not found\n", num_cell);
        }

        neighbors[face] = id_neighbor_cell * CELL_SIZE + id_neighbor_face;
        neighbors[id_neighbor_cell * CELL_SIZE + id_neighbor_face] = face;
      } else if (idc->GetNumberOfIds() == 0)
        neighbors[face] = e_neigh_code_out_bound; // граничная ячейка
      else {
        RETURN_ERR("More than 1 neighbor????\n");
      }
    }
  }

  return e_completion_success;
}

static void NormalAndSquareFace2D(int cell_number, int face_number, const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, Type &S, Vector3 &n) {

  vtkIdList *idp = unstructured_grid->GetCell(cell_number)->GetEdge(face_number)->GetPointIds();

  double P0[3], P1[3], P3[3]; // узлы ячейки

  unstructured_grid->GetPoint(idp->GetId(0), P0);
  unstructured_grid->GetPoint(idp->GetId(1), P1);

  double a[2];
  for (int i = 0; i < 2; i++) {
    a[i] = P1[i] - P0[i];
  }
  n[0] = P0[1] - P1[1];
  n[1] = P1[0] - P0[0];
  n[2] = 0;

  S = sqrt(a[0] * a[0] + a[1] * a[1]); // длина отрезка

  // нормировка
  n[0] /= S;
  n[1] /= S;

  vtkSmartPointer<vtkIdList> idp2 = unstructured_grid->GetCell(cell_number)->GetPointIds();

  int id; // номер узла не ячейки не принадлежащей рассматриваемой гране
  for (int i = 0; i < CELL_SIZE; i++) {
    int count = 0;
    for (int j = 0; j < 2; j++)
      if (idp2->GetId(i) != idp->GetId(j))
        count++;
    if (count == 2) {
      id = i;
      break;
    }
  }

  //    определение ориентации нормали
  Type sum = 0;
  unstructured_grid->GetPoint(idp2->GetId(id), P3);
  for (int i = 0; i < 2; i++) {
    sum += n[i] * (P3[i] - P0[i]);
  }

  if (sum < 0)
    for (int i = 0; i < 2; i++)
      n[i] *= -1;
}

int GetNormalAndAreas2D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Normals> &normals, std::vector<Type> &areas) {

  const int n = unstructured_grid->GetNumberOfCells();
  normals.resize(n, Normals(CELL_SIZE));
  areas.resize(n * (CELL_SIZE));

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < CELL_SIZE; j++) {
      NormalAndSquareFace2D(i, j, unstructured_grid, areas[i * CELL_SIZE + j], normals[i].n[j]);
    }
  }

  return e_completion_success;
}

int GetVolume2D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Type> &volumes) {

  auto GetCellVolume{[](size_t cell_number, const vtkSmartPointer<vtkUnstructuredGrid> &u_grid) {
    double P0[3], P1[3], P2[3];

    vtkSmartPointer<vtkIdList> idp = u_grid->GetCell(cell_number)->GetPointIds();
    u_grid->GetPoint(idp->GetId(0), P0);
    u_grid->GetPoint(idp->GetId(1), P1);
    u_grid->GetPoint(idp->GetId(2), P2);

    Eigen::Vector2d a, b;
    for (int i = 0; i < 2; i++) {
      a[i] = P1[i] - P0[i];
      b[i] = P2[i] - P0[i];
    }
    return (Eigen::Vector3d(a[0], a[1], 0).cross(Eigen::Vector3d(b[0], b[1], 0))).norm() / 2;
  }};

  const int n = unstructured_grid->GetNumberOfCells();
  volumes.resize(n);

  for (int i = 0; i < n; i++) {
    volumes[i] = GetCellVolume(i, unstructured_grid);
  }

  return e_completion_success;
}

static void CenterOfTetra2D(const int number_cell, const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, Vector3 &point_in_tetra) {

  double P0[3], P1[3], P2[3];
  vtkSmartPointer<vtkIdList> idp = unstructured_grid->GetCell(number_cell)->GetPointIds();

  unstructured_grid->GetPoint(idp->GetId(0), P0);
  unstructured_grid->GetPoint(idp->GetId(1), P1);
  unstructured_grid->GetPoint(idp->GetId(2), P2);

  for (int k = 0; k < 2; k++) {
    point_in_tetra[k] = (P0[k] + P1[k] + P2[k]) / 3;
  }
  point_in_tetra[2] = 0;
}

int GetCentersOfCells2D(const vtkSmartPointer<vtkUnstructuredGrid> &unstructured_grid, std::vector<Vector3> &centers) {

  const int n = unstructured_grid->GetNumberOfCells();
  centers.resize(n);

  for (int i = 0; i < n; i++) {
    CenterOfTetra2D(i, unstructured_grid, centers[i]);
  }

  return e_completion_success;
}

#endif // NUMBER_OF_MEASUREMENTS

#endif //! USE_VTK
