#ifndef VECTOR_H
#define VECTOR_H

#include "type.hpp"
#include <cmath>

template <typename T>
class   Vector {
public:
    T   x, y, z;

    Vector(T _x = 0, T _y = 0, T _z = 0) : x(_x), y(_y), z(_z) {}
    Vector(const Vector<T> &vec) : x(vec.x), y(vec.y), z(vec.z) {}

    T   len() const { return sqrt(x*x + y*y + z*z); }
    T   len2() const { return x*x + y*y + z*z; }
    T   sum() const { return x + y + z; }

    Vector<T>&  operator= (const Vector<T> &vec);
    Vector<T>&  operator+=(const Vector<T> &vec);
    Vector<T>&  operator-=(const Vector<T> &vec);
    Vector<T>&  operator^=(const Vector<T> &vec);
    Vector<T>&  operator+=(const T &val);
    Vector<T>&  operator-=(const T &val);
    Vector<T>&  operator*=(const T &val);
    Vector<T>&  operator/=(const T &val);

//    friend Vector<T>    operator+(const Vector<T> &vec);
//    friend Vector<T>    operator-(const Vector<T> &vec);

//    // Vector and Vector
//    friend Vector<T>    operator+(const Vector<T> &vec1, const Vector<T> &vec2);
//    friend Vector<T>    operator-(const Vector<T> &vec1, const Vector<T> &vec2);
//    friend Vector<T>    operator^(const Vector<T> &vec1, const Vector<T> &vec2);    // cross product
//    friend T            operator*(const Vector<T> &vec1, const Vector<T> &vec2);    // dot product

//    // Vector and Value
//    friend Vector<T>    operator+(const Vector<T> &vec, const T &val);
//    friend Vector<T>    operator-(const Vector<T> &vec, const T &val);
//    friend Vector<T>    operator*(const Vector<T> &vec, const T &val);
//    friend Vector<T>    operator/(const Vector<T> &vec, const T &val);
//    friend Vector<T>    operator+(const T &val, const Vector<T> &vec);
//    friend Vector<T>    operator-(const T &val, const Vector<T> &vec);
//    friend Vector<T>    operator*(const T &val, const Vector<T> &vec);
//    friend Vector<T>    operator/(const T &val, const Vector<T> &vec);

    // overload iostream operator
    // TODO
};


template <typename T>
Vector<T>&  Vector<T>::operator=(const Vector<T> &vec)
{
    x = vec.x;
    y = vec.y;
    z = vec.z;
    return *this;
}

template <typename T>
Vector<T>&  Vector<T>::operator+=(const Vector<T> &vec)
{
    (*this) = (*this) + vec;
    return *this;
}

template <typename T>
Vector<T>&  Vector<T>::operator-=(const Vector<T> &vec)
{
    (*this) = (*this) - vec;
    return *this;
}

template <typename T>
Vector<T>&  Vector<T>::operator^=(const Vector<T> &vec)
{
    (*this) = (*this) ^ vec;
    return *this;
}

template <typename T>
Vector<T>&  Vector<T>::operator+=(const T &val)
{
    (*this) = (*this) + val;
    return *this;
}

template <typename T>
Vector<T>&  Vector<T>::operator-=(const T &val)
{
    (*this) = (*this) - val;
    return *this;
}

template <typename T>
Vector<T>&  Vector<T>::operator*=(const T &val)
{
    (*this) = (*this) * val;
    return *this;
}

template <typename T>
Vector<T>&  Vector<T>::operator/=(const T &val)
{
    (*this) = (*this) / val;
    return *this;
}

template <typename T>
Vector<T>   operator+(const Vector<T> &vec)
{ return  Vector<T>(vec); }

template <typename T>
Vector<T>   operator-(const Vector<T> &vec)
{ return  Vector<T>(-vec.x, -vec.y, -vec.z); }

template <typename T>
Vector<T>   operator+(const Vector<T> &vec1, const Vector<T> &vec2)
{ return Vector<T>(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z); }

template <typename T>
Vector<T>   operator-(const Vector<T> &vec1, const Vector<T> &vec2)
{ return Vector<T>(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z); }

template <typename T>
Vector<T>   operator^(const Vector<T> &vec1, const Vector<T> &vec2)
{
    return Vector<T>( vec1.y * vec2.z - vec1.z * vec2.y,
                      vec1.z * vec2.x - vec1.x * vec2.z,
                      vec1.x * vec2.y - vec1.y * vec2.x );
}

template <typename T>
T   operator*(const Vector<T> &vec1, const Vector<T> &vec2)
{ return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z; }

template <typename T>
Vector<T>   operator+(const Vector<T> &vec, const T &val)
{ return Vector<T>(vec.x + val, vec.y + val, vec.z + val); }

template <typename T>
Vector<T>   operator-(const Vector<T> &vec, const T &val)
{ return Vector<T>(vec.x - val, vec.y - val, vec.z - val); }

template <typename T>
Vector<T>   operator*(const Vector<T> &vec, const T &val)
{ return Vector<T>(vec.x * val, vec.y * val, vec.z * val); }

template <typename T>
Vector<T>   operator/(const Vector<T> &vec, const T &val)
{ return Vector<T>(vec.x / val, vec.y / val, vec.z / val); }

template <typename T>
Vector<T>   operator+(const T &val, const Vector<T> &vec)
{ return vec + val; }

template <typename T>
Vector<T>   operator-(const T &val, const Vector<T> &vec)
{ return (vec - val) * -1; }

template <typename T>
Vector<T>   operator*(const T &val, const Vector<T> &vec)
{ return vec * val; }

template <typename T>
Vector<T>   operator/(const T &val, const Vector<T> &vec)
{ return Vector<T>(val / vec.x, val / vec.y, val / vec.z); }


typedef  Vector<INTG>   VecIntg;
typedef  Vector<REAL>   VecReal;

#endif // VECTOR_H
