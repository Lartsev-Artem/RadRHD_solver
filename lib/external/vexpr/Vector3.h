#ifndef VECTOR3_H
#define VECTOR3_H

namespace sf
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

        static constexpr Vector3<T> Zero() {return Vector3<T>(0,0,0);}
        
        ////////////////////////////////////////////////////////////
        /// \brief Default constructor
        ///
        /// Creates a `Vector3(0, 0, 0)`.
        ///
        ////////////////////////////////////////////////////////////
        constexpr Vector3() :_data{0,0,0}{}// = default;
        //constexpr Vector3()  = default;

        ////////////////////////////////////////////////////////////
        /// \brief Construct the vector from its coordinates
        ///
        /// \param x X coordinate
        /// \param y Y coordinate
        /// \param z Z coordinate
        ///
        ////////////////////////////////////////////////////////////
        constexpr Vector3(T x, T y, T z);
        constexpr Vector3(const Vector3& x): _data{x.x(),x.y(),x.z()}{}
        constexpr Vector3(const T* ptr): _data{ptr[0],ptr[1],ptr[2]}{}

        ////////////////////////////////////////////////////////////
        /// \brief Converts the vector to another type of vector
        ///
        ////////////////////////////////////////////////////////////
        template <typename U>
        constexpr explicit operator Vector3<U>() const;

        
        constexpr T operator[](const int i) const
        {
            //assert
            return _data[i];            
        }
        constexpr T& operator[](const int i)
        {
            //assert
            return _data[i];
        }

        const T* data() const {return _data;}

        ////////////////////////////////////////////////////////////
        /// \brief Length of the vector <i><b>(floating-point)</b></i>.
        ///
        /// If you are not interested in the actual length, but only in comparisons, consider using `lengthSquared()`.
        ///
        ////////////////////////////////////////////////////////////
        [[nodiscard]] T length() const;

        ////////////////////////////////////////////////////////////
        /// \brief Square of vector's length.
        ///
        /// Suitable for comparisons, more efficient than `length()`.
        ///
        ////////////////////////////////////////////////////////////
        [[nodiscard]] constexpr T lengthSquared() const;

        ////////////////////////////////////////////////////////////
        /// \brief Vector with same direction but length 1 <i><b>(floating-point)</b></i>.
        ///
        /// \pre `*this` is no zero vector.
        ///
        ////////////////////////////////////////////////////////////
        [[nodiscard]] Vector3 normalized() const;

        ////////////////////////////////////////////////////////////
        /// \brief Dot product of two 3D vectors.
        ///
        ////////////////////////////////////////////////////////////
        [[nodiscard]] constexpr T dot(const Vector3 &rhs) const;

        ////////////////////////////////////////////////////////////
        /// \brief Cross product of two 3D vectors.
        ///
        ////////////////////////////////////////////////////////////
        [[nodiscard]] constexpr Vector3 cross(const Vector3 &rhs) const;

        ////////////////////////////////////////////////////////////
        /// \brief Component-wise multiplication of `*this` and `rhs`.
        ///
        /// Computes `(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z)`.
        ///
        /// Scaling is the most common use case for component-wise multiplication/division.
        /// This operation is also known as the Hadamard or Schur product.
        ///
        ////////////////////////////////////////////////////////////
        [[nodiscard]] constexpr Vector3 componentWiseMul(const Vector3 &rhs) const;

        ////////////////////////////////////////////////////////////
        /// \brief Component-wise division of `*this` and `rhs`.
        ///
        /// Computes `(lhs.x/rhs.x, lhs.y/rhs.y, lhs.z/rhs.z)`.
        ///
        /// Scaling is the most common use case for component-wise multiplication/division.
        ///
        /// \pre Neither component of `rhs` is zero.
        ///
        ////////////////////////////////////////////////////////////
        [[nodiscard]] constexpr Vector3 componentWiseDiv(const Vector3 &rhs) const;

        ////////////////////////////////////////////////////////////
        // Member data
        ////////////////////////////////////////////////////////////        
        constexpr T x() const {return _data[0];} //!< X coordinate of the vector
        constexpr T y() const {return _data[1];} //!< Y coordinate of the vector
        constexpr T z() const {return _data[2];} //!< Z coordinate of the vector        

        T _data[3]; //!< XYZ coordinate of the vector
    };

    ////////////////////////////////////////////////////////////
    /// \relates Vector3
    /// \brief Overload of unary `operator-`
    ///
    /// \param left Vector to negate
    ///
    /// \return Member-wise opposite of the vector
    ///
    ////////////////////////////////////////////////////////////
    template <typename T>
    [[nodiscard]] constexpr Vector3<T> operator-(const Vector3<T> &left);

    ////////////////////////////////////////////////////////////
    /// \relates Vector3
    /// \brief Overload of binary `operator+=`
    ///
    /// This operator performs a member-wise addition of both vectors,
    /// and assigns the result to `left`.
    ///
    /// \param left  Left operand (a vector)
    /// \param right Right operand (a vector)
    ///
    /// \return Reference to `left`
    ///
    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> &operator+=(Vector3<T> &left, const Vector3<T> &right);

    ////////////////////////////////////////////////////////////
    /// \relates Vector3
    /// \brief Overload of binary `operator-=`
    ///
    /// This operator performs a member-wise subtraction of both vectors,
    /// and assigns the result to `left`.
    ///
    /// \param left  Left operand (a vector)
    /// \param right Right operand (a vector)
    ///
    /// \return Reference to `left`
    ///
    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> &operator-=(Vector3<T> &left, const Vector3<T> &right);

    ////////////////////////////////////////////////////////////
    /// \relates Vector3
    /// \brief Overload of binary `operator+`
    ///
    /// \param left  Left operand (a vector)
    /// \param right Right operand (a vector)
    ///
    /// \return Member-wise addition of both vectors
    ///
    ////////////////////////////////////////////////////////////
    template <typename T>
    [[nodiscard]] constexpr Vector3<T> operator+(const Vector3<T> &left, const Vector3<T> &right);

    ////////////////////////////////////////////////////////////
    /// \relates Vector3
    /// \brief Overload of binary `operator-`
    ///
    /// \param left  Left operand (a vector)
    /// \param right Right operand (a vector)
    ///
    /// \return Member-wise subtraction of both vectors
    ///
    ////////////////////////////////////////////////////////////
    template <typename T>
    [[nodiscard]] constexpr Vector3<T> operator-(const Vector3<T> &left, const Vector3<T> &right);

    ////////////////////////////////////////////////////////////
    /// \relates Vector3
    /// \brief Overload of binary `operator*`
    ///
    /// \param left  Left operand (a vector)
    /// \param right Right operand (a scalar value)
    ///
    /// \return Member-wise multiplication by `right`
    ///
    ////////////////////////////////////////////////////////////
    template <typename T>
    [[nodiscard]] constexpr Vector3<T> operator*(const Vector3<T> &left, T right);

    ////////////////////////////////////////////////////////////
    /// \relates Vector3
    /// \brief Overload of binary `operator*`
    ///
    /// \param left  Left operand (a scalar value)
    /// \param right Right operand (a vector)
    ///
    /// \return Member-wise multiplication by `left`
    ///
    ////////////////////////////////////////////////////////////
    template <typename T>
    [[nodiscard]] constexpr Vector3<T> operator*(T left, const Vector3<T> &right);

    ////////////////////////////////////////////////////////////
    /// \relates Vector3
    /// \brief Overload of binary `operator*=`
    ///
    /// This operator performs a member-wise multiplication by `right`,
    /// and assigns the result to `left`.
    ///
    /// \param left  Left operand (a vector)
    /// \param right Right operand (a scalar value)
    ///
    /// \return Reference to `left`
    ///
    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> &operator*=(Vector3<T> &left, T right);

    ////////////////////////////////////////////////////////////
    /// \relates Vector3
    /// \brief Overload of binary `operator/`
    ///
    /// \param left  Left operand (a vector)
    /// \param right Right operand (a scalar value)
    ///
    /// \return Member-wise division by `right`
    ///
    ////////////////////////////////////////////////////////////
    template <typename T>
    [[nodiscard]] constexpr Vector3<T> operator/(const Vector3<T> &left, T right);

    ////////////////////////////////////////////////////////////
    /// \relates Vector3
    /// \brief Overload of binary `operator/=`
    ///
    /// This operator performs a member-wise division by `right`,
    /// and assigns the result to `left`.
    ///
    /// \param left  Left operand (a vector)
    /// \param right Right operand (a scalar value)
    ///
    /// \return Reference to `left`
    ///
    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> &operator/=(Vector3<T> &left, T right);

    ////////////////////////////////////////////////////////////
    /// \relates Vector3
    /// \brief Overload of binary `operator==`
    ///
    /// This operator compares strict equality between two vectors.
    ///
    /// \param left  Left operand (a vector)
    /// \param right Right operand (a vector)
    ///
    /// \return `true` if `left` is equal to `right`
    ///
    ////////////////////////////////////////////////////////////
    template <typename T>
    [[nodiscard]] constexpr bool operator==(const Vector3<T> &left, const Vector3<T> &right);

    ////////////////////////////////////////////////////////////
    /// \relates Vector3
    /// \brief Overload of binary `operator!=`
    ///
    /// This operator compares strict difference between two vectors.
    ///
    /// \param left  Left operand (a vector)
    /// \param right Right operand (a vector)
    ///
    /// \return `true` if `left` is not equal to `right`
    ///
    ////////////////////////////////////////////////////////////
    template <typename T>
    [[nodiscard]] constexpr bool operator!=(const Vector3<T> &left, const Vector3<T> &right);

    // Aliases for the most common types
    using Vector3i = Vector3<int>;
    using Vector3f = Vector3<float>;
    using Vector3d = Vector3<double>;

} // namespace sf

#include "Vector3.inl"

#endif //! VECTOR3_H