
#include "Vector3.h"

#include <assert.h>

namespace sf
{
////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T>::Vector3(T x, T y, T z) : _data{x,y,z}
    {
    }
    ////////////////////////////////////////////////////////////
    template <typename T>
    template <typename U>
    constexpr Vector3<T>::operator Vector3<U>() const
    {
        return Vector3<U>(static_cast<U>(x()), static_cast<U>(y()), static_cast<U>(z()));
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr T Vector3<T>::lengthSquared() const
    {
        return dot(*this);
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr T Vector3<T>::dot(const Vector3<T> &rhs) const
    {
        return x() * rhs.x() + y() * rhs.y() + z() * rhs.z();
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> Vector3<T>::cross(const Vector3<T> &rhs) const
    {
        return Vector3<T>((y()* rhs.z()) - (z() * rhs.y()), (z() * rhs.x()) - (x() * rhs.z()), (x() * rhs.y()) - (y()* rhs.x()));
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> Vector3<T>::componentWiseMul(const Vector3<T> &rhs) const
    {
        return Vector3<T>(x() * rhs.x(), y()* rhs.y(), z() * rhs.z());
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> Vector3<T>::componentWiseDiv(const Vector3<T> &rhs) const
    {
        assert(rhs.x() != 0 && "Vector3::componentWiseDiv() cannot divide by 0");
        assert(rhs.y() != 0 && "Vector3::componentWiseDiv() cannot divide by 0");
        assert(rhs.z() != 0 && "Vector3::componentWiseDiv() cannot divide by 0");
        return Vector3<T>(x() / rhs.x(), y() / rhs.y(), z() / rhs.z());
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> operator-(const Vector3<T> &left)
    {
        return Vector3<T>(-left.x(), -left.y(), -left.z());
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> &operator+=(Vector3<T> &left, const Vector3<T> &right)
    {
        left.x() += right.x();
        left.y() += right.y();
        left.z() += right.z();

        return left;
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> &operator-=(Vector3<T> &left, const Vector3<T> &right)
    {
        left.x() -= right.x();
        left.y() -= right.y();
        left.z() -= right.z();

        return left;
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> operator+(const Vector3<T> &left, const Vector3<T> &right)
    {
        return Vector3<T>(left.x() + right.x(), left.y() + right.y(), left.z() + right.z());
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> operator-(const Vector3<T> &left, const Vector3<T> &right)
    {
        return Vector3<T>(left.x() - right.x(), left.y() - right.y(), left.z() - right.z());
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> operator*(const Vector3<T> &left, T right)
    {
        return Vector3<T>(left.x() * right, left.y() * right, left.z() * right);
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> operator*(T left, const Vector3<T> &right)
    {
        return Vector3<T>(right.x() * left, right.y() * left, right.z() * left);
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> &operator*=(Vector3<T> &left, T right)
    {
        left.x() *= right;
        left.y() *= right;
        left.z() *= right;

        return left;
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> operator/(const Vector3<T> &left, T right)
    {
        assert(right != 0 && "Vector3::operator/ cannot divide by 0");
        return Vector3<T>(left.x() / right, left.y() / right, left.z() / right);
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr Vector3<T> &operator/=(Vector3<T> &left, T right)
    {
        assert(right != 0 && "Vector3::operator/= cannot divide by 0");
        left.x() /= right;
        left.y() /= right;
        left.z() /= right;

        return left;
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr bool operator==(const Vector3<T> &left, const Vector3<T> &right)
    {
        return (left.x() == right.x()) && (left.y() == right.y()) && (left.z() == right.z());
    }

    ////////////////////////////////////////////////////////////
    template <typename T>
    constexpr bool operator!=(const Vector3<T> &left, const Vector3<T> &right)
    {
        return !(left == right);
    }

} // namespace sf