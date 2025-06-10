#include "../Matrix3.h"

#include <assert.h>

namespace vexpr {
namespace internal {
template <typename T>
constexpr TransposedMatrix3<T>::TransposedMatrix3(const Matrix3<T> &A) : _ref_matrix(A) {
}

template <typename T>
constexpr T TransposedMatrix3<T>::operator()(size_t i, size_t j) const noexcept {
  return _ref_matrix(j, i);
}
}; // namespace internal

//!-------------------------- Constructs --------------------------//
template <typename T>
constexpr Matrix3<T>::Matrix3() noexcept
    : _data{0, 0, 0, 0, 0, 0, 0, 0, 0} {
}

template <typename T>
constexpr Matrix3<T>::Matrix3(const Matrix3 &v) noexcept
    : _data{v.data()} {
}

template <typename T>
constexpr Matrix3<T>::Matrix3(const T *p) {
  std::copy(p, p + 9, _data);
}

template <typename T>
constexpr Matrix3<T>::Matrix3(std::initializer_list<T> init) {
  std::copy(init.begin(), init.end(), _data);
}
//!-------------------------- Methods --------------------------//
template <typename T>
constexpr const T *Matrix3<T>::data() const noexcept {
  return _data;
}

template <typename T>
std::string Matrix3<T>::serialize() const {
  std::string str = "Matrix3:";
  for (size_t i = 0; i < 3; i++) {
    str += "\t";
    for (size_t j = 0; j < 2; j++) {
      str += std::to_string(_data[3 * i + j]) + ", ";
    }
    str += std::to_string(_data[3 * i + 2]) + "\n\t";
  }
  return str + "\n";
}

template <typename T>
constexpr internal::TransposedMatrix3<T> Matrix3<T>::transpose() const {
  return internal::TransposedMatrix3<T>(*this);
}

//!-------------------------- Operators --------------------------//
template <typename T>
template <typename U>
constexpr Matrix3<T>::operator Matrix3<U>() const {
  return Matrix3<U>({
      static_cast<U>(_data[0]),
      static_cast<U>(_data[1]),
      static_cast<U>(_data[2]),
      static_cast<U>(_data[3]),
      static_cast<U>(_data[4]),
      static_cast<U>(_data[5]),
      static_cast<U>(_data[6]),
      static_cast<U>(_data[7]),
      static_cast<U>(_data[8]),
  });
}

template <typename T>
constexpr T Matrix3<T>::operator()(const size_t i, const size_t j) const noexcept {
  assert(i < 3 && "Matrix3::operator[] const. i > =3");
  assert(i < 3 && "Matrix3::operator[] const. j > =3");
  return _data[3 * i + j];
}

template <typename T>
constexpr T &Matrix3<T>::operator()(const size_t i, const size_t j) noexcept {
  assert(i < 3 && "Matrix3::operator(). i > =3");
  assert(i < 3 && "Matrix3::operator(). j > =3");
  return _data[3 * i + j];
}

template <typename T>
constexpr const T *Matrix3<T>::operator[](const size_t i) const noexcept {
  assert(i < 3 && "Matrix3::operator[] const idx > =3");
  return _data + 3 * i;
}

template <typename T>
constexpr T *Matrix3<T>::operator[](const size_t i) noexcept {
  assert(i < 3 && "Matrix3::operator[] idx > =3");
  return _data + 3 * i;
}

////////////////////////////////////////////////////////////

template <typename T, typename MatrixExpr>
constexpr Vector3<T> operator*(const MatrixExpr &A, const Vector3<T> &b) noexcept {
  return Vector3<T>{
      A(0, 0) * b[0] + A(0, 1) * b[1] + A(0, 2) * b[2],
      A(1, 0) * b[0] + A(1, 1) * b[1] + A(1, 2) * b[2],
      A(2, 0) * b[0] + A(2, 1) * b[1] + A(2, 2) * b[2]};
}

////////////////////////////////////////////////////////////
template <typename T>
constexpr bool operator==(const Matrix3<T> &left, const Matrix3<T> &right) noexcept {
  bool res = true;
  for (size_t i = 0; i < 9; i++) {
    res = (res && (left._data[i] == right._data[i]));
  }
  return res;
}

template <typename T>
constexpr bool operator!=(const Matrix3<T> &left, const Matrix3<T> &right) noexcept {
  return !(left == right);
}

} // namespace vexpr