#if defined BUILD_GRAPH && defined USE_CUDA
#include "graph_inner_bound.h"

#include "global_value.h"

#include "reader_txt.h"
#include "writer_bin.h"

#include "cuda_ray_interface.h"
#include "linear_alg.h"

int graph::trace_through_boundary::InitDevice()
{

  std::vector<FaceCell> inter_faces; ///< внутренние грани

  if (files_sys::txt::ReadSimple(glb_files.base_address + F_FACE_ID, inter_faces))
    return e_completion_fail;

  std::vector<Ray_t> rays(inter_faces.size());
  cuda::ray_tracing::interface::InitDevice(inter_faces, rays);

  return e_completion_success;
}

void graph::trace_through_boundary::ClearDevice()
{
  cuda::ray_tracing::interface::ClearDevice();
}

int graph::trace_through_boundary::FindBoundCondOnInnerBoundary(int num_dir, const Vector3 &direction, const std::map<IntId, FaceCell> &faces, const std::set<IntId> &outer_part,
                                                                std::vector<IntId> &intersections)
{

  std::vector<Ray_t> rays(outer_part.size());
  intersections.resize(outer_part.size());

  int k = 0;
  for (auto out : outer_part)
  {
    const FaceCell f = faces.find(out)->second;
    Vector3 center = (f.face.A + f.face.B + f.face.C) / 3.;
    rays[k].orig = center - direction * (1e-4); // сдвиг чтобы не было самопересечения
    rays[k].direction = -direction;             // ищем обратное пересечение
    k++;
  }

  return cuda::ray_tracing::interface::FindInnerIntersection(rays, intersections);
}

void graph::trace_through_boundary::SortInnerBoundary(const std::vector<IntId> &graph, IntersectBound_t &intersections)
{

  size_t N = intersections.out_id_cell.size();

  std::vector<std::pair<int, int>> pairs(N);
  for (size_t i = 0; i < N; i++)
  {
    pairs[i] = std::make_pair(intersections.out_id_cell[i], intersections.code[i]);
  }

  // переупорядочивание определяющих граней на внутренней границе в соответствии с их порядком в решении
  std::vector<int> buf(N, 0);
  int pos = 0;
  for (auto num_cell : graph)
  {
    for (size_t j = 0; j < N; j++)
    {
      if (num_cell == pairs[j].first)
      {
        intersections.code[pos++] = pairs[j].second;
        break;
      }
    }
  }

#if 0 // def RRHD_DEBUG
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
      DIE_IF(intersections.out_id_cell[i] == intersections.code[j] / 4); //ячейка не должна определяться сама на себя
    }

#endif
}

#endif //! MAKE_TRACE