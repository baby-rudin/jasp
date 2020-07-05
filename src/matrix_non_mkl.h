#ifndef MATRIX_NON_MKL_H
#define MATRIX_NON_MKL_H

#include "type.hpp"
#include "constant.h"
#include <iostream>
#include <string>

class Matrix {
private:
    size_t      _row;
    size_t      _col;
    REAL        *_data;

public:
    Matrix(size_t row = 0, size_t col = 1, REAL val = 0.0);
    Matrix(size_t row, size_t col, REAL *array);
    Matrix(size_t row, size_t col, REAL **array);
    Matrix(const Matrix &mat);
    ~Matrix(void);

    size_t  rows() const { return _row; }
    size_t  cols() const { return _col; }
    size_t  size() const { return _row * _col; }

    std::string     print() const;

    Matrix  inverse() const;
    Matrix  trans() const;

    REAL&   operator()(size_t pos)
    { return _data[pos]; }

    REAL  operator()(size_t pos) const
    { return _data[pos]; }

    REAL&   operator()(size_t row, size_t col)
    { return _data[row * _col + col]; }

    REAL  operator()(size_t row, size_t col) const
    { return _data[row * _col + col]; }

    Matrix& operator= (const Matrix &mat);
    Matrix& operator+=(const Matrix &mat);
    Matrix& operator-=(const Matrix &mat);
    Matrix& operator*=(const Matrix &mat);
    Matrix& operator+=(REAL val);
    Matrix& operator-=(REAL val);
    Matrix& operator*=(REAL val);
    Matrix& operator/=(REAL val);

    friend Matrix   operator+(const Matrix &mat);
    friend Matrix   operator-(const Matrix &mat);

    // Matrix and Matrix
    friend Matrix   operator+(const Matrix &mat1, const Matrix &mat2);
    friend Matrix   operator-(const Matrix &mat1, const Matrix &mat2);
    friend Matrix   operator*(const Matrix &mat1, const Matrix &mat2);

    // Matrix and value
    friend Matrix   operator+(const Matrix &mat, REAL val);
    friend Matrix   operator-(const Matrix &mat, REAL val);
    friend Matrix   operator*(const Matrix &mat, REAL val);
    friend Matrix   operator/(const Matrix &mat, REAL val);
    friend Matrix   operator+(REAL val, const Matrix &mat);
    friend Matrix   operator-(REAL val, const Matrix &mat);
    friend Matrix   operator*(REAL val, const Matrix &mat);
    friend Matrix   operator/(REAL val, const Matrix &mat);

    // overload ostream operator
    friend std::ostream &operator<<(std::ostream &os, const Matrix &mat);

    friend class EigenSolver;
};


class EigenSolver {
private:
    size_t  _dim;
    Matrix  _EigenVals;
    Matrix  _EigenVecs;

public:
    EigenSolver(const Matrix &mat);

    Matrix  eig_val() const { return _EigenVals; }
    Matrix  eig_vec() const { return _EigenVecs; }
};


#endif // MATRIX_H
