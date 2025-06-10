/**
 * @file VectorX.h
 * @author Artem
 * @brief constexpr VectorN class
 * @version 0.1
 * @date 2025-06-10
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef VECTORX_H
#define VECTORX_H

#include <cstddef>
#include <initializer_list>
#include <string>

namespace vexpr
{
  ////////////////////////////////////////////////////////////
  /// \brief Utility template class for manipulating
  ///        3-dimensional vectors
  ///
  ////////////////////////////////////////////////////////////
  template <typename T, size_t _N = 3>
  class VectorX
  {
  public:
    static constexpr VectorX<T, _N> Zero() noexcept { return VectorX<T, _N>(); }

    //!-------------------------- Constructs --------------------------//

    /// \brief Construct the default vector
    constexpr VectorX() noexcept;

    /// \brief Construct the vector from its coordinates
    constexpr VectorX(std::initializer_list<T> init);

    /// \brief Construct the vector from other vector
    constexpr VectorX(const VectorX &v) noexcept;

    /// \brief Construct the vector from array data
    /// \param p - array of length >= 3
    constexpr VectorX(const T *p);

    //!-------------------------- Methods --------------------------//
    [[nodiscard]] constexpr const T *data() const noexcept;

    /// \brief Norm of the vector
    /// \warning Not constexpr method!!!
    [[nodiscard]] T norm() const;

    /// \brief Vector with same direction but length 1
    /// \warning Not constexpr method!!!
    [[nodiscard]] VectorX normalized() const;

    ///\brief Serialize data of vector to string
    /// \warning Not constexpr method!!!
    std::string serialize() const;

    /// \brief Norm^2 of the vector
    /// \note Better to use, then norm()
    [[nodiscard]] constexpr T norm2() const noexcept;

    /// \brief Dot product of two 3D vectors.
    [[nodiscard]] constexpr T dot(const VectorX &rhs) const noexcept;

    //!-------------------------- Operators --------------------------//

    /// \brief Converts the vector to another type of vector
    template <typename U, size_t N2>
    constexpr explicit operator VectorX<U, N2>() const;

    constexpr T operator[](const size_t i) const noexcept;
    constexpr T &operator[](const size_t i) noexcept;

    //!------------------------- Member data -------------------------//
  private:
    T _data[_N]; ///< XYZ coordinate of the vector
  };

  ////////////////////////////////////////////////////////////
  /// \relates VectorX
  /// \brief Overload of unary `operator-`
  template <typename T>
  [[nodiscard]] constexpr VectorX<T> operator-(const VectorX<T> &left) noexcept;

  ////////////////////////////////////////////////////////////
  /// \relates VectorX
  /// \brief Overload of binary `operator+`
  template <typename T>
  [[nodiscard]] constexpr VectorX<T> operator+(const VectorX<T> &left, const VectorX<T> &right) noexcept;

  /// \relates VectorX
  /// \brief Overload of binary `operator-`
  template <typename T>
  [[nodiscard]] constexpr VectorX<T> operator-(const VectorX<T> &left, const VectorX<T> &right) noexcept;

  /// \relates VectorX
  /// \brief Overload of binary `operator+`
  template <typename T>
  [[nodiscard]] constexpr VectorX<T> operator/(const VectorX<T> &left, const VectorX<T> &right) noexcept;

  /// \relates VectorX
  /// \brief Overload of binary `operator-`
  template <typename T>
  [[nodiscard]] constexpr VectorX<T> operator*(const VectorX<T> &left, const VectorX<T> &right) noexcept;

  ////////////////////////////////////////////////////////////

  /// \relates VectorX
  /// \brief Overload of binary `operator*`
  template <typename T>
  [[nodiscard]] constexpr VectorX<T> operator*(const VectorX<T> &left, T val) noexcept;

  /// \relates VectorX
  /// \brief Overload of binary `operator/`
  template <typename T>
  [[nodiscard]] constexpr VectorX<T> operator/(const VectorX<T> &left, T val);

  /// \relates VectorX
  /// \brief Overload of binary `operator*`
  template <typename T>
  [[nodiscard]] constexpr VectorX<T> operator*(T val, const VectorX<T> &right) noexcept;

  ////////////////////////////////////////////////////////////

  /// \relates VectorX
  template <typename T>
  constexpr void operator*=(VectorX<T> &left, T val) noexcept;

  /// \relates VectorX
  /// \brief Overload of binary `operator/=`
  template <typename T>
  constexpr void operator/=(VectorX<T> &left, T val);

  /// \relates VectorX
  /// \brief Overload of binary `operator+=`
  template <typename T>
  constexpr void operator+=(VectorX<T> &left, const VectorX<T> &right) noexcept;

  /// \relates VectorX
  /// \brief Overload of binary `operator-=`
  template <typename T>
  constexpr void operator-=(VectorX<T> &left, const VectorX<T> &right) noexcept;

  ////////////////////////////////////////////////////////////

  /// \relates VectorX
  /// \brief Overload of binary `operator==`
  template <typename T>
  [[nodiscard]] constexpr bool operator==(const VectorX<T> &left, const VectorX<T> &right) noexcept;

  /// \relates VectorX
  /// \brief Overload of binary `operator!=`
  template <typename T>
  [[nodiscard]] constexpr bool operator!=(const VectorX<T> &left, const VectorX<T> &right) noexcept;

  ////////////////////////////////////////////////////////////

} // namespace vexpr

#include "inl/VectorX.inl"

#endif //! VectorX_H