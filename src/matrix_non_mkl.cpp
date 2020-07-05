#include "matrix_non_mkl.h"
#include "eigen_sym.h"
#include <iostream>
#include <string>
#include <algorithm>


#define BUFF_LEN 1024

using namespace std;

// ==================== Matrix ====================

Matrix::Matrix(size_t row, size_t col, REAL val)
{
    _row    = row;
    _col    = col;

    if (size())
        _data = new REAL [size()];
    else
        _data = nullptr;

    for (size_t i=0; i<size(); i++)
        _data[i] = val;
}

Matrix::Matrix(size_t row, size_t col, REAL *array)
{
    _row    = row;
    _col    = col;

    if (size())
        _data = new REAL [size()];
    else
        _data = nullptr;

    for (size_t i=0; i<size(); i++)
        _data[i] = array[i];
}

Matrix::Matrix(size_t row, size_t col, REAL **array)
{
    _row    = row;
    _col    = col;

    if (size())
        _data = new REAL [size()];
    else
        _data = nullptr;

    for (size_t i=0; i<_row; i++)
        for (size_t j=0; j<_col; j++)
            (*this)(i,j) = array[i][j];
}

Matrix::Matrix(const Matrix &mat)
{
    _row    = mat._row;
    _col    = mat._col;

    if (size())
        _data = new REAL [size()];
    else
        _data = nullptr;

    for (size_t i=0; i<size(); i++)
        _data[i] = mat._data[i];
}

Matrix::~Matrix(void)
{
    delete [] _data;
    _data = nullptr;
}


string  Matrix::print() const
{
    string  ret;
    char    buff[BUFF_LEN];

    for (size_t i=0; i<_row; i++) {
        for (size_t j=0; j<_col; j++) {
            sprintf(buff, "%14.6lf", (*this)(i, j));
            ret += string(buff);
        }
        ret += "\n";
    }

    return ret;
}


Matrix  Matrix::inverse() const
{
    if (_row != _col) {
        cout << "calculate inverse for row != col matrix" << endl;
        exit(234);
    }

    Matrix ret(*this);

    int *ipiv = new int [_row];
//    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, int(_row), int(_row), ret._data, int(_row), ipiv);
//    LAPACKE_dgetri(LAPACK_ROW_MAJOR, int(_row), ret._data, int(_row), ipiv);
    delete [] ipiv;
    ipiv = nullptr;

    return  ret;
}

Matrix  Matrix::trans() const
{
    Matrix ret(_col, _row);
    for (size_t i = 0; i < ret._row; i++)
        for (size_t j = 0; j < ret._col; j++)
            ret(i, j) = (*this)(j, i);

    return ret;
}

// operators of Matrix
Matrix  &Matrix::operator= (const Matrix &mat)
{
    if (this == &mat) return *this;

    _row    = mat._row;
    _col    = mat._col;
    delete [] _data;
    _data   = new REAL [size()];
    for (size_t i=0; i<size(); i++)
        _data[i] = mat._data[i];

    return *this;
}

Matrix  &Matrix::operator+=(const Matrix &mat)
{
    (*this) = (*this) + mat;
    return *this;
}

Matrix  &Matrix::operator-=(const Matrix &mat)
{
    (*this) = (*this) - mat;
    return *this;
}

Matrix  &Matrix::operator*=(const Matrix &mat)
{
    (*this) = (*this) * mat;
    return *this;
}

Matrix  &Matrix::operator+=(REAL val)
{
    for (size_t i=0; i<size(); i++)
        _data[i] += val;
    return *this;
}

Matrix  &Matrix::operator-=(REAL val)
{
    for (size_t i=0; i<size(); i++)
        _data[i] -= val;
    return *this;
}

Matrix  &Matrix::operator*=(REAL val)
{
    for (size_t i=0; i<size(); i++)
        _data[i] *= val;
    return *this;
}

Matrix  &Matrix::operator/=(REAL val)
{
    for (size_t i=0; i<size(); i++)
        _data[i] /= val;
    return *this;
}

Matrix  operator+(const Matrix &mat)
{ return Matrix(mat); }

Matrix  operator-(const Matrix &mat)
{
    Matrix  ret(mat);

    for (size_t i=0; i<ret.size(); i++)
        ret._data[i] = -ret._data[i];

    return  ret;
}

// operator: Matrix and Matrix
Matrix  operator+(const Matrix &mat1, const Matrix &mat2)
{
    if (mat1._row != mat2._row || mat1._col != mat2._col) {
        cout << "+- two matrix with different size" << endl;
        exit(201);
    }

    Matrix  ret(mat1);
    for (size_t i=0; i<mat1.size(); i++)
        ret._data[i] += mat2._data[i];

    return  ret;
}

Matrix  operator-(const Matrix &mat1, const Matrix &mat2)
{
    if (mat1._row != mat2._row || mat1._col != mat2._col) {
        cout << "+- two matrix with different size" << endl;
        exit(201);
    }

    Matrix  ret(mat1);
    for (size_t i=0; i<mat1.size(); i++)
        ret._data[i] -= mat2._data[i];

    return  ret;
}

Matrix  operator*(const Matrix &mat1, const Matrix &mat2)
{
    if (mat1._col != mat2._row) {
        cout << "* two matrix, while mat1.col != mat2.row" << endl;
    }

    Matrix  ret(mat1._row, mat2._col);
    for (size_t i=0; i<ret._row; i++)
        for (size_t j=0; j<ret._col; j++)
            for (size_t k=0; k<mat1._col; k++)
                ret(i,j) += mat1(i,k) * mat2(k,j);

    return  ret;
}

// operator: Matrix and value
Matrix  operator+(const Matrix &mat, REAL val)
{
    Matrix  ret(mat);
    for (size_t i=0; i<mat.size(); i++)
        ret._data[i] += val;

    return ret;
}

Matrix  operator-(const Matrix &mat, REAL val)
{
    Matrix  ret(mat);
    for (size_t i=0; i<mat.size(); i++)
        ret._data[i] -= val;

    return ret;
}

Matrix  operator*(const Matrix &mat, REAL val)
{
    Matrix  ret(mat);
    for (size_t i=0; i<mat.size(); i++)
        ret._data[i] *= val;

    return ret;
}

Matrix  operator/(const Matrix &mat, REAL val)
{
    Matrix  ret(mat);
    for (size_t i=0; i<mat.size(); i++)
        ret._data[i] /= val;

    return ret;
}

Matrix  operator+(REAL val, const Matrix &mat)
{ return  mat + val; }

Matrix  operator-(REAL val, const Matrix &mat)
{ return (mat - val) * REAL(-1.0); }

Matrix  operator*(REAL val, const Matrix &mat)
{ return mat * val; }

Matrix  operator/(REAL val, const Matrix &mat)
{
    Matrix  ret(mat);

    for (size_t i=0; i<mat.size(); i++)
        ret._data[i] = val / mat._data[i];

    return ret;
}

ostream &operator<<(ostream &os, const Matrix &mat)
{
    os << mat.print() << endl;
    return os;
}

// =================== EigenSolver ====================

EigenSolver::EigenSolver(const Matrix &mat)
{
    if (mat._row != mat._col) {
        cout << "calculate eigen for row != col matrix" << endl;
        exit(203);
    }

    _dim        = mat._row;
    _EigenVals  = Matrix(_dim, _dim);
    _EigenVecs  = Matrix(_dim, _dim);

    REAL  *wr     = new REAL [_dim];
    REAL  *vr     = new REAL [mat.size()];
    REAL  *src    = new REAL [mat.size()];
    for (size_t i=0; i<mat.size(); i++)
        src[i] = mat._data[i];

    if (!eig_sym(src, INTG(_dim), wr, vr)) {
        cout << "eigensolver didn't converge!" << endl;
        exit(1353);
    }

    // sort by eigenvalues
    size_t *id = new size_t [_dim];
    for (size_t i = 0; i < _dim; i++)
        id[i] = i;
    sort(id, id + _dim, [wr](const int a, const int b) {
        return wr[a] < wr[b];
    });

    for (size_t i = 0; i < _dim; i++) {
        _EigenVals(i, i) = wr[id[i]];
        for (size_t j = 0; j < _dim; j++)
            _EigenVecs(j, i) = vr[ j * _dim + id[i] ];
    }

    delete [] wr;   wr  = nullptr;
    delete [] vr;   vr  = nullptr;
    delete [] src;  src = nullptr;
    delete [] id;   id  = nullptr;
}
