#ifndef MATHFUN_H
#define MATHFUN_H

#include <string>
#include <iostream>

class Vec;          // class Vec used to represent a vector
class Mat;          // class Mat used to represent a matrix
class EigenSolver;  // calculate eigenvalues and eigenvectors of Mat


class Vec
{
private:
    double  x, y, z;

public:
    Vec(double x = 0.0);
    Vec(double x, double y, double z);

    double          len() const;
    double          len2() const;
    std::string     print() const;

    Vec    &operator= (const Vec &vector);
    Vec    &operator+=(const Vec &vector);
    Vec    &operator-=(const Vec &vector);
    Vec    &operator+=(double Val);
    Vec    &operator-=(double Val);
    Vec    &operator*=(double Val);
    Vec    &operator/=(double Val);

public:
    friend Vec      operator+(const Vec &vector);
    friend Vec      operator-(const Vec &vector);
    friend Vec      operator+(const Vec &vector1, const Vec &vector2);
    friend Vec      operator-(const Vec &vector1, const Vec &vector2);
    friend Vec      operator^(const Vec &vector1, const Vec &vector2);  // cross product
    friend double   operator*(const Vec &vector1, const Vec &vector2);  // dot product

    friend Vec      operator+(const Vec &vector, double Val);
    friend Vec      operator-(const Vec &vector, double Val);
    friend Vec      operator*(const Vec &vector, double Val);
    friend Vec      operator/(const Vec &vector, double Val);
    friend Vec      operator+(double Val, const Vec &vector);
    friend Vec      operator-(double Val, const Vec &vector);
    friend Vec      operator*(double Val, const Vec &vector);
    friend Vec      operator/(double Val, const Vec &vector);

    friend std::istream &operator>>(std::istream &is, Vec &vector);
    friend std::ostream &operator<<(std::ostream &os, const Vec &vector);
};


class Mat
{
private:
    int     Row;
    int     Col;
    double  *Data;      // row major array

public:
    Mat(void);
    Mat(int Rows, int Cols, double Val = 0.0);
    Mat(int Rows, int Cols, double *Array);
    Mat(int Rows, int Cols, double **Array);
    Mat(const Mat &matrix);

    ~Mat(void);

    void            reshape(int Rows, int Cols);
    std::string     print() const;

    int     rows() const
    {
        return Row;
    }
    int     cols() const
    {
        return Col;
    }
    int     size() const
    {
        return Row * Col;
    }

    // matrix function
    Mat     inv();
    Mat     trans();

    double &operator()(int R, int C);
    double  operator()(int R, int C) const;   // which const ???

    Mat    &operator= (const Mat &matrix);
    Mat    &operator+=(const Mat &matrix);
    Mat    &operator-=(const Mat &matrix);
    Mat    &operator*=(const Mat &matrix);
    Mat    &operator+=(double Val);
    Mat    &operator-=(double Val);
    Mat    &operator*=(double Val);

public:
    friend Mat  operator+(const Mat &matrix);
    friend Mat  operator-(const Mat &matrix);
    friend Mat  operator+(const Mat &matrix1, const Mat &matrix2);
    friend Mat  operator-(const Mat &matrix1, const Mat &matrix2);
    friend Mat  operator*(const Mat &matrix1, const Mat &matrix2);
    friend Mat  operator+(const Mat &matrix, double Val);
    friend Mat  operator-(const Mat &matrix, double Val);
    friend Mat  operator*(const Mat &matrix, double Val);
    friend Mat  operator+(double Val, const Mat &matrix);
    friend Mat  operator-(double Val, const Mat &matrix);
    friend Mat  operator*(double Val, const Mat &matrix);

    friend std::ostream &operator<<(std::ostream &os, const Mat &matrix);

    friend class EigenSolver;
};


class EigenSolver
{
private:
    int     Dim;
    Mat     EigenVals;
    Mat     EigenVecs;

public:
    EigenSolver(const Mat &matrix);

    Mat     eig_val() const
    {
        return EigenVals;
    }
    Mat     eig_vec() const
    {
        return EigenVecs;
    }

};


long long   factor(int n);    // calculate the factor

#endif // MATHFUN_H
