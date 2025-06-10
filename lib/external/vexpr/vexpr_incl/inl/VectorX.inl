#include "../VectorX.h"

#include <assert.h>

namespace vexpr {
//!-------------------------- Constructs --------------------------//

template <typename T, size_t N>
constexpr VectorX<T, N>::VectorX() noexcept
    : _data{} {
}

template <typename T, size_t N>
constexpr VectorX<T, N>::VectorX(std::initializer_list<T> init)
    : _data{} {
  std::copy(init.begin(), init.end(), _data);
}

template <typename T, size_t N>
constexpr VectorX<T, N>::VectorX(const VectorX &v) noexcept
    : _data{} {
  std::copy(v._data, v._data + N, _data);
}

template <typename T, size_t N>
constexpr VectorX<T, N>::VectorX(const T *p)
    : _data{} {
  std::copy(p, p + N, _data);
}
//!-------------------------- Methods --------------------------//

template <typename T, size_t N>
constexpr const T *VectorX<T, N>::data() const noexcept {
  return _data;
}
////////////////////////////////////////////////////////////
template <typename T, size_t N>
T VectorX<T, N>::norm() const {
  return sqrt(this->dot(*this));
}
template <typename T, size_t N>
VectorX<T, N> VectorX<T, N>::normalized() const {
  assert(norm2() > 1e-300 && "VectorX::normalized() zero vector!!\n");
  return *this / norm();
}

template <typename T, size_t N>
std::string VectorX<T, N>::serialize() const {
  std::string str = "VectorX:\t";
  for (size_t i = 0; i < N - 1; i++) {
    str += std::to_string(_data[i]) + ", ";
  }
  return str + std::to_string(_data[N - 1]) + "\n";
}

template <typename T, size_t N>
constexpr T VectorX<T, N>::norm2() const noexcept {
  return this->dot(*this);
}

template <typename T, size_t N>
constexpr T VectorX<T, N>::dot(const VectorX<T, N> &rhs) const noexcept {

  T sum = 0;
  for (size_t i = 0; i < N; i++) {
    sum[i] += this->_data[i] * rhs[i];
  }
  return sum;
}

////////////////////////////////////////////////////////////

//!-------------------------- Operators --------------------------//
template <typename T, size_t N>
constexpr T VectorX<T, N>::operator[](const size_t i) const noexcept {
  assert(i < N && "VectorX::operator[] const idx > =N");
  return _data[i];
}

template <typename T, size_t N>
constexpr T &VectorX<T, N>::operator[](const size_t i) noexcept {
  assert(i < N && "VectorX::operator[] idx > =N");
  return _data[i];
}

////////////////////////////////////////////////////////////

template <typename T, size_t N>
constexpr VectorX<T, N> operator-(const VectorX<T, N> &left) noexcept {

  VectorX<T, N> x;
  for (size_t i = 0; i < N; i++) {
    x[i] = -left[i];
  }
  return x;
}

////////////////////////////////////////////////////////////

template <typename T, size_t N>
constexpr VectorX<T, N> operator+(const VectorX<T, N> &left, const VectorX<T, N> &right) noexcept {

  VectorX<T, N> x;
  for (size_t i = 0; i < N; i++) {
    x[i] = left[i] + right[i];
  }
  return x;
}

template <typename T, size_t N>
constexpr VectorX<T, N> operator-(const VectorX<T, N> &left, const VectorX<T, N> &right) noexcept {
  VectorX<T, N> x;
  for (size_t i = 0; i < N; i++) {
    x[i] = left[i] - right[i];
  }
  return x;
}

template <typename T, size_t N>
constexpr VectorX<T, N> operator/(const VectorX<T, N> &left, const VectorX<T, N> &right) noexcept {

  VectorX<T, N> x;
  for (size_t i = 0; i < N; i++) {
    assert(right[i] != 0 && "VectorX::componentWiseDiv() cannot divide by 0");
    x[i] = left[i] / right[i];
  }
  return x;
}

template <typename T, size_t N>
constexpr VectorX<T, N> operator*(const VectorX<T, N> &left, const VectorX<T, N> &right) noexcept {
  VectorX<T, N> x;
  for (size_t i = 0; i < N; i++) {
    x[i] = left[i] * right[i];
  }
  return x;
}

////////////////////////////////////////////////////////////
template <typename T, size_t N>
constexpr VectorX<T, N> operator*(const VectorX<T, N> &left, T val) noexcept {
  VectorX<T, N> x;
  for (size_t i = 0; i < N; i++) {
    x[i] = left[i] * val;
  }
  return x;
}

template <typename T, size_t N>
constexpr VectorX<T, N> operator*(T val, const VectorX<T, N> &right) noexcept {
  VectorX<T, N> x;
  for (size_t i = 0; i < N; i++) {
    x[i] = val * right[i];
  }
  return x;
}

template <typename T, size_t N>
constexpr VectorX<T, N> operator/(const VectorX<T, N> &left, T val) {
  assert(val != 0 && "VectorX::operator/ cannot divide by 0");
  VectorX<T, N> x;
  for (size_t i = 0; i < N; i++) {
    x[i] = left[i] / val;
  }
  return x;
}

////////////////////////////////////////////////////////////

template <typename T, size_t N>
constexpr void operator*=(VectorX<T, N> &left, T val) noexcept {
  for (size_t i = 0; i < N; i++) {
    x[i] = left[i] *= val;
  }
}

template <typename T, size_t N>
constexpr void operator/=(VectorX<T, N> &left, T val) {
  assert(val != 0 && "VectorX::operator/= cannot divide by 0");
  for (size_t i = 0; i < N; i++) {
    x[i] = left[i] /= val;
  }
}

template <typename T, size_t N>
constexpr void operator+=(VectorX<T, N> &left, const VectorX<T, N> &right) noexcept {
  for (size_t i = 0; i < N; i++) {
    left[i] += right[i];
  }
}

template <typename T, size_t N>
constexpr void operator-=(VectorX<T, N> &left, const VectorX<T, N> &right) noexcept {
  for (size_t i = 0; i < N; i++) {
    left[i] -= right[i];
  }
}

////////////////////////////////////////////////////////////
template <typename T, size_t N>
constexpr bool operator==(const VectorX<T, N> &left, const VectorX<T, N> &right) noexcept {
  bool res = true;
  for (size_t i = 0; i < 9; i++) {
    res = (res && (left[i] == right[i]));
  }
  return res;
}

template <typename T, size_t N>
constexpr bool operator!=(const VectorX<T, N> &left, const VectorX<T, N> &right) noexcept {
  return !(left == right);
}

} // namespace vexpr