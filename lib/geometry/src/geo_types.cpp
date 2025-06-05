
#include "geo_types.h"

std::ostream &operator<<(std::ostream &os, const FaceCell &f)
{
  return os << std::setprecision(16)
            << f.face_id << ' '
            << f.face.A[0] << ' ' << f.face.A[1] << ' ' << f.face.A[2] << ' '
            << f.face.B[0] << ' ' << f.face.B[1] << ' ' << f.face.B[2] << ' '
            << f.face.C[0] << ' ' << f.face.C[1] << ' ' << f.face.C[2];
}
std::ostream &operator<<(std::ostream &os, const std::pair<const int, FaceCell> &f)
{
  return os << f.second;
}

std::istream &operator>>(std::istream &is, FaceCell &f)
{
  return is >> f.face_id >>
         f.face.A[0] >> f.face.A[1] >> f.face.A[2] >>
         f.face.B[0] >> f.face.B[1] >> f.face.B[2] >>
         f.face.C[0] >> f.face.C[1] >> f.face.C[2];
}