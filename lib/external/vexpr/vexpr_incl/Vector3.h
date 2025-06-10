/**
 * @file Vector3.h
 * @author Artem
 * @brief constexpr Vector3 class
 * @version 0.1
 * @date 2025-06-10
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef VECTOR3_H
#define VECTOR3_H

#include <cstddef>
#include <string>

namespace vexpr
{
    ////////////////////////////////////////////////////////////
    /// \brief Utility template class for manipulating
    ///        3-dimensional vectors
    ///
    ////////////////////////////////////////////////////////////
    template <typename T>
    class Vector3
    {
    public:
        static constexpr Vector3<T> Zero() noexcept { return Vector3<T>(); }

        //!-------------------------- Constructs --------------------------//

        /// \brief Construct the vector from its coordinates
        constexpr Vector3(T x = 0, T y = 0, T z = 0) noexcept;

        /// \brief Construct the vector from other vector
        constexpr Vector3(const Vector3 &v) noexcept;

        /// \brief Construct the vector from array data
        /// \param p - array of length >= 3
        constexpr Vector3(const T *p);

        //!-------------------------- Methods --------------------------//

        [[nodiscard]] constexpr T x() const noexcept;
        [[nodiscard]] constexpr T y() const noexcept;
        [[nodiscard]] constexpr T z() const noexcept;
        [[nodiscard]] constexpr const T *data() const noexcept;

        /// \brief Norm of the vector
        /// \warning Not constexpr method!!!
        [[nodiscard]] T norm() const;

        /// \brief Vector with same direction but length 1
        /// \warning Not constexpr method!!!
        [[nodiscard]] Vector3 normalized() const;

        ///\brief Serialize data of vector to string
        /// \warning Not constexpr method!!!
        std::string serialize() const;

        /// \brief Norm^2 of the vector
        /// \note Better to use, then norm()
        [[nodiscard]] constexpr T norm2() const noexcept;

        /// \brief Dot product of two 3D vectors.
        [[nodiscard]] constexpr T dot(const Vector3 &rhs) const noexcept;

        /// \brief Cross product of two 3D vectors.
        [[nodiscard]] constexpr Vector3 cross(const Vector3 &rhs) const noexcept;

        //!-------------------------- Operators --------------------------//

        /// \brief Converts the vector to another type of vector
        template <typename U>
        constexpr explicit operator Vector3<U>() const;

        constexpr T operator[](const size_t i) const noexcept;
        constexpr T &operator[](const size_t i) noexcept;

        //!------------------------- Member data -------------------------//
    private:
        T _data[3]; ///< XYZ coordinate of the vector
    };

    ////////////////////////////////////////////////////////////
    /// \relates Vector3
    /// \brief Overload of unary `operator-`
    template <typename T>
    [[nodiscard]] constexpr Vector3<T> operator-(const Vector3<T> &left) noexcept;

    ////////////////////////////////////////////////////////////
    /// \relates Vector3
    /// \brief Overload of binary `operator+`
    template <typename T>
    [[nodiscard]] constexpr Vector3<T> operator+(const Vector3<T> &left, const Vector3<T> &right) noexcept;

    /// \relates Vector3
    /// \brief Overload of binary `operator-`
    template <typename T>
    [[nodiscard]] constexpr Vector3<T> operator-(const Vector3<T> &left, const Vector3<T> &right) noexcept;

    /// \relates Vector3
    /// \brief Overload of binary `operator+`
    template <typename T>
    [[nodiscard]] constexpr Vector3<T> operator/(const Vector3<T> &left, const Vector3<T> &right) noexcept;

    /// \relates Vector3
    /// \brief Overload of binary `operator-`
    template <typename T>
    [[nodiscard]] constexpr Vector3<T> operator*(const Vector3<T> &left, const Vector3<T> &right) noexcept;

    ////////////////////////////////////////////////////////////

    /// \relates Vector3
    /// \brief Overload of binary `operator*`
    template <typename T>
    [[nodiscard]] constexpr Vector3<T> operator*(const Vector3<T> &left, T val) noexcept;

    /// \relates Vector3
    /// \brief Overload of binary `operator/`
    template <typename T>
    [[nodiscard]] constexpr Vector3<T> operator/(const Vector3<T> &left, T val);

    /// \relates Vector3
    /// \brief Overload of binary `operator*`
    template <typename T>
    [[nodiscard]] constexpr Vector3<T> operator*(T val, const Vector3<T> &right) noexcept;

    ////////////////////////////////////////////////////////////

    /// \relates Vector3
    template <typename T>
    constexpr void operator*=(Vector3<T> &left, T val) noexcept;

    /// \relates Vector3
    /// \brief Overload of binary `operator/=`
    template <typename T>
    constexpr void operator/=(Vector3<T> &left, T val);

    /// \relates Vector3
    /// \brief Overload of binary `operator+=`
    template <typename T>
    constexpr void operator+=(Vector3<T> &left, const Vector3<T> &right) noexcept;

    /// \relates Vector3
    /// \brief Overload of binary `operator-=`
    template <typename T>
    constexpr void operator-=(Vector3<T> &left, const Vector3<T> &right) noexcept;

    ////////////////////////////////////////////////////////////

    /// \relates Vector3
    /// \brief Overload of binary `operator==`
    template <typename T>
    [[nodiscard]] constexpr bool operator==(const Vector3<T> &left, const Vector3<T> &right) noexcept;

    /// \relates Vector3
    /// \brief Overload of binary `operator!=`
    template <typename T>
    [[nodiscard]] constexpr bool operator!=(const Vector3<T> &left, const Vector3<T> &right) noexcept;

    ////////////////////////////////////////////////////////////

    // Aliases for the most common types
    using Vector3i = Vector3<int>;
    using Vector3f = Vector3<float>;
    using Vector3d = Vector3<double>;

} // namespace vexpr

#include "inl/Vector3.inl"

#endif //! VECTOR3_H