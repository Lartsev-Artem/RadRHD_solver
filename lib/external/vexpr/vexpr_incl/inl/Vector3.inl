#include "../Vector3.h"

#include <assert.h>

namespace vexpr
{
    //!-------------------------- Constructs --------------------------//
    template <typename T>
    constexpr Vector3<T>::Vector3(T x, T y, T z) noexcept
        : _data{x, y, z}
    {
    }

    template <typename T>
    constexpr Vector3<T>::Vector3(const Vector3 &v) noexcept
        : _data{v.x(), v.y(), v.z()}
    {
    }

    template <typename T>
    constexpr Vector3<T>::Vector3(const T *p)
        : _data{p[0], p[1], p[2]}
    {
    }
    //!-------------------------- Methods --------------------------//
    template <typename T>
    constexpr T Vector3<T>::x() const noexcept
    {
        return _data[0];
    }

    template <typename T>
    constexpr T Vector3<T>::y() const noexcept
    {
        return _data[1];
    }

    template <typename T>
    constexpr T Vector3<T>::z() const noexcept
    {
        return _data[2];
    }

    template <typename T>
    constexpr const T *Vector3<T>::data() const noexcept
    {
        return _data;
    }
    ////////////////////////////////////////////////////////////
    template <typename T>
    T Vector3<T>::norm() const
    {
        return sqrt(this->dot(*this));
    }
    template <typename T>
    Vector3<T> Vector3<T>::normalized() const
    {
        assert(norm2() > 1e-300 && "Vector3::normalized() zero vector!!\n");
        return *this / norm();
    }

    template <typename T>
    std::string Vector3<T>::serialize() const
    {
        std::string str = "Vector3:\t";
        for (size_t i = 0; i < 2; i++)
        {
            str += std::to_string(_data[i]) + ", ";
        }
        return str + std::to_string(_data[2]) + "\n";
    }

    template <typename T>
    constexpr T Vector3<T>::norm2() const noexcept
    {
        return this->dot(*this);
    }

    template <typename T>
    constexpr T Vector3<T>::dot(const Vector3<T> &rhs) const noexcept
    {
        return x() * rhs.x() + y() * rhs.y() + z() * rhs.z();
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> Vector3<T>::cross(const Vector3<T> &rhs) const noexcept
    {
        return Vector3<T>((y() * rhs.z()) - (z() * rhs.y()),
                          (z() * rhs.x()) - (x() * rhs.z()),
                          (x() * rhs.y()) - (y() * rhs.x()));
    }

    //!-------------------------- Operators --------------------------//
    template <typename T>
    template <typename U>
    constexpr Vector3<T>::operator Vector3<U>() const
    {
        return Vector3<U>(static_cast<U>(x()), static_cast<U>(y()), static_cast<U>(z()));
    }

    template <typename T>
    constexpr T Vector3<T>::operator[](const size_t i) const noexcept
    {
        assert(i < 3 && "Vector3::operator[] const idx > =3");
        return _data[i];
    }

    template <typename T>
    constexpr T &Vector3<T>::operator[](const size_t i) noexcept
    {
        assert(i < 3 && "Vector3::operator[] idx > =3");
        return _data[i];
    }

    ////////////////////////////////////////////////////////////

    template <typename T>
    constexpr Vector3<T> operator-(const Vector3<T> &left) noexcept
    {
        return Vector3<T>(-left.x(), -left.y(), -left.z());
    }

    ////////////////////////////////////////////////////////////

    template <typename T>
    constexpr Vector3<T> operator+(const Vector3<T> &left, const Vector3<T> &right) noexcept
    {
        return Vector3<T>(left.x() + right.x(),
                          left.y() + right.y(),
                          left.z() + right.z());
    }

    template <typename T>
    constexpr Vector3<T> operator-(const Vector3<T> &left, const Vector3<T> &right) noexcept
    {
        return Vector3<T>(left.x() - right.x(),
                          left.y() - right.y(),
                          left.z() - right.z());
    }

    template <typename T>
    constexpr Vector3<T> operator/(const Vector3<T> &left, const Vector3<T> &right) noexcept
    {
        assert(right.x() != 0 && "Vector3::componentWiseDiv() cannot divide by 0");
        assert(right.y() != 0 && "Vector3::componentWiseDiv() cannot divide by 0");
        assert(right.z() != 0 && "Vector3::componentWiseDiv() cannot divide by 0");
        return Vector3<T>(left.x() / right.x(),
                          left.y() / right.y(),
                          left.z() / right.z());
    }

    template <typename T>
    constexpr Vector3<T> operator*(const Vector3<T> &left, const Vector3<T> &right) noexcept
    {
        return Vector3<T>(left.x() * right.x(),
                          left.y() * right.y(),
                          left.z() * right.z());
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> operator*(const Vector3<T> &left, T val) noexcept
    {
        return Vector3<T>(left.x() * val,
                          left.y() * val,
                          left.z() * val);
    }

    template <typename T>
    constexpr Vector3<T> operator*(T val, const Vector3<T> &right) noexcept
    {
        return Vector3<T>(val * right.x(),
                          val * right.y(),
                          val * right.z());
    }

    template <typename T>
    constexpr Vector3<T> operator/(const Vector3<T> &left, T val)
    {
        assert(val != 0 && "Vector3::operator/ cannot divide by 0");
        return Vector3<T>(left.x() / val,
                          left.y() / val,
                          left.z() / val);
    }

    ////////////////////////////////////////////////////////////

    template <typename T>
    constexpr void operator*=(Vector3<T> &left, T val) noexcept
    {
        left.x() *= val;
        left.y() *= val;
        left.z() *= val;
    }

    template <typename T>
    constexpr void operator/=(Vector3<T> &left, T val)
    {
        assert(val != 0 && "Vector3::operator/= cannot divide by 0");
        left.x() /= val;
        left.y() /= val;
        left.z() /= val;
    }

    template <typename T>
    constexpr void operator+=(Vector3<T> &left, const Vector3<T> &right) noexcept
    {
        left.x() += right.x();
        left.y() += right.y();
        left.z() += right.z();
    }

    template <typename T>
    constexpr void operator-=(Vector3<T> &left, const Vector3<T> &right) noexcept
    {
        left.x() -= right.x();
        left.y() -= right.y();
        left.z() -= right.z();
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr bool operator==(const Vector3<T> &left, const Vector3<T> &right) noexcept
    {
        return (left.x() == right.x()) &&
               (left.y() == right.y()) &&
               (left.z() == right.z());
    }

    template <typename T>
    constexpr bool operator!=(const Vector3<T> &left, const Vector3<T> &right) noexcept
    {
        return !(left == right);
    }

} // namespace vexpr