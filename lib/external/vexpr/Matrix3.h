#ifndef MATRIX3_H
#define MATRIX3_H

#include "Vector3.h"
#include <array>

namespace sf
{
    ////////////////////////////////////////////////////////////
    /// \brief Utility template class for manipulating
    ///        3-dimensional matrixes
    ///
    ////////////////////////////////////////////////////////////
template <typename T, int _size = 3>  class Matrix3;

template <typename T, int _size = 3>
class TransposedMatrix3 {
public:
    constexpr TransposedMatrix3() = delete;
    constexpr TransposedMatrix3(const TransposedMatrix3& m) = delete;
    constexpr explicit TransposedMatrix3(const Matrix3<T,_size>& mat) : matrix(mat) {}    

    constexpr T operator()(int i, int j) const {  return matrix(j, i);  }
    
private:
    const Matrix3<T,_size>& matrix;
};

    template <typename T, int _size>
    class Matrix3
    {
    public:
        constexpr const T* data() const {return _data;}        

        //constexpr Matrix3() = default;
        // constexpr Matrix3() :_data{}{            
        //        for (int i = 0; i < 3; ++i) _data[i] = 0;
        
        // }

        // // Инициализация списком (row-major)
        constexpr Matrix3(std::initializer_list<T> init) {
        std::copy(init.begin(), init.end(), _data);
        }

        static constexpr Matrix3<T> Zero() {return Matrix3();}
    
        // Возврат прокси-объекта для транспонирования (без копирования)
        constexpr auto transpose() const { return TransposedMatrix3<T,_size>(*this); }

        constexpr T operator()(const int i,const int j) const
        {
            //assert
            return _data[_size*i+j];
        }
        constexpr T& operator()(const int i, const int j)
        {
            //assert
            return _data[_size*i+j];
        }

        constexpr const T* operator[](const int i) const
        {
            //assert
            return _data + _size*i;
        }
        constexpr T* operator[](const int i)
        {
            //assert
            return _data + _size*i;
        }

        ////////////////////////////////////////////////////////////
        // Member data
        ////////////////////////////////////////////////////////////
        T _data[_size *_size];
    };
    
// Умножение любой матричного "выражения" на Vector3
template <typename T, typename MatrixExpr>
constexpr Vector3<T> operator*(const MatrixExpr& mat, const Vector3<T>& vec) {
        return Vector3<T>{
        mat(0,0)*vec[0] + mat(0,1)*vec[1] + mat(0,2)*vec[2],
        mat(1,0)*vec[0] + mat(1,1)*vec[1] + mat(1,2)*vec[2],
        mat(2,0)*vec[0] + mat(2,1)*vec[1] + mat(2,2)*vec[2]
        };
}

   
    // Aliases for the most common types
    using Matrix3i = Matrix3<int,3>;
    using Matrix3f = Matrix3<float,3>;
    using Matrix3d = Matrix3<double,3>;

} // namespace sf

//#include <Matrix3.inl>

#endif //! Matrix3_H