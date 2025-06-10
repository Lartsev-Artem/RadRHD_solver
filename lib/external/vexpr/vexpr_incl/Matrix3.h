/**
 * @file Matrix3.h
 * @author Artem
 * @brief constexpr Matrix3 class
 * @version 0.1
 * @date 2025-06-10
 *
 * @copyright Copyright (c) 2025
 *
 */
#ifndef MATRIX3_H
#define MATRIX3_H

#include "Vector3.h"
#include <initializer_list>

namespace vexpr
{
    template <typename T>
    class Matrix3;

    namespace internal
    {
        /// \brief the transposed representation
        template <typename T>
        class TransposedMatrix3
        {
        public:
            constexpr TransposedMatrix3() = delete;
            constexpr TransposedMatrix3(const TransposedMatrix3 &A) = delete;
            constexpr explicit TransposedMatrix3(const Matrix3<T> &A);

            constexpr T operator()(size_t i, size_t j) const noexcept;

        private:
            const Matrix3<T> &_ref_matrix;
        };
    };

    ////////////////////////////////////////////////////////////
    /// \brief Utility template class for manipulating
    ///        3-dimensional matrixes
    ///
    ////////////////////////////////////////////////////////////
    template <typename T>
    class Matrix3
    {
    public:
        static constexpr Matrix3<T> Zero() noexcept { return Matrix3<T>(); }

        //!-------------------------- Constructs -------------------------//
        /// \brief Construct the default matrix
        constexpr Matrix3() noexcept;

        /// \brief Construct the matrix from other matrix
        constexpr Matrix3(const Matrix3 &v) noexcept;

        /// \brief Construct the vector from array data
        /// \param p - array of length >= 9
        constexpr Matrix3(const T *p);

        ///\brief Initialization with a list (row-major)
        constexpr Matrix3(std::initializer_list<T> init);

        //!--------------------------- Methods ---------------------------//
        constexpr const T *data() const noexcept;

        ///\brief Serialize data of vector to string
        /// \warning Not constexpr method!!!
        std::string serialize() const;

        ///\brief Return of the proxy object for transposition (without copying)
        constexpr internal::TransposedMatrix3<T> transpose() const;

        //!-------------------------- Operators --------------------------//

        /// \brief Converts the vector to another type of vector
        template <typename U>
        constexpr explicit operator Matrix3<U>() const;

        constexpr T operator()(const size_t i, const size_t j) const noexcept;
        constexpr T &operator()(const size_t i, const size_t j) noexcept;
        constexpr const T *operator[](const size_t i) const noexcept;
        constexpr T *operator[](const size_t i) noexcept;

        //!------------------------- Member data -------------------------//
    private:
        T _data[9];
    };

    /// \relates Matrix3
    /// \brief Overload of binary `operator==`
    template <typename T>
    [[nodiscard]] constexpr bool operator==(const Matrix3<T> &left, const Matrix3<T> &right) noexcept;

    /// \relates Matrix3
    /// \brief Overload of binary `operator!=`
    template <typename T>
    [[nodiscard]] constexpr bool operator!=(const Matrix3<T> &left, const Matrix3<T> &right) noexcept;

    /**
     * \relates Matrix3
     * \brief Overload of binary `operator*`
     * \details Multiplication of any matrix "expression" by Vector3
     */
    template <typename T, typename MatrixExpr>
    constexpr Vector3<T> operator*(const MatrixExpr &A, const Vector3<T> &b) noexcept;

    // Aliases for the most common types
    using Matrix3i = Matrix3<int>;
    using Matrix3f = Matrix3<float>;
    using Matrix3d = Matrix3<double>;

} // namespace vexpr

#include "inl/Matrix3.inl"

#endif //! Matrix3_H