#include "mathfun.h"
#include <string>
#include <cmath>
#include <algorithm>

#include "mkl.h"


#define BUFF_LEN 1024

using namespace std;

//==================== Vec ====================

Vec::Vec(double x)
{
    this->x = x;
    this->y = x;
    this->z = x;
}

Vec::Vec(double x, double y, double z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

double  Vec::len() const
{
    return sqrt((*this) * (*this));
}

double  Vec::len2() const
{
    return (*this) * (*this);
}

string  Vec::print() const
{
    char buff[BUFF_LEN];
    sprintf(buff, "%14.8lf%14.8lf%14.8lf",
            x, y, z);
    return  string(buff);
}


// operator of Vec
Vec    &Vec::operator= (const Vec &vector)
{
    if (this == &vector)
        return *this;

    x = vector.x;
    y = vector.y;
    z = vector.z;

    return *this;
}

Vec    &Vec::operator+=(const Vec &vector)
{
    x += vector.x;
    y += vector.y;
    z += vector.z;
    return *this;
}

Vec    &Vec::operator-=(const Vec &vector)
{
    x -= vector.x;
    y -= vector.y;
    z -= vector.z;
    return *this;
}

Vec    &Vec::operator+=(double Val)
{
    x += Val;
    y += Val;
    z += Val;
    return *this;
}

Vec    &Vec::operator-=(double Val)
{
    x -= Val;
    y -= Val;
    z -= Val;
    return *this;
}

Vec    &Vec::operator*=(double Val)
{
    x *= Val;
    y *= Val;
    z *= Val;
    return *this;
}

Vec    &Vec::operator/=(double Val)
{
    x /= Val;
    y /= Val;
    z /= Val;
    return *this;
}


// friend operator of Vec

Vec     operator+(const Vec &vector)
{
    return Vec(vector);
}

Vec     operator-(const Vec &vector)
{
    return Vec(-vector.x,
               -vector.y,
               -vector.z);
}

Vec     operator+(const Vec &vector1, const Vec &vector2)
{
    return Vec(vector1.x + vector2.x,
               vector1.y + vector2.y,
               vector1.z + vector2.z);
}

Vec     operator-(const Vec &vector1, const Vec &vector2)
{
    return Vec(vector1.x - vector2.x,
               vector1.y - vector2.y,
               vector1.z - vector2.z);
}

Vec     operator^(const Vec &V1, const Vec &V2)     // cross product
{
    return Vec(V1.y * V2.z - V1.z * V2.y,
               V1.z * V2.x - V1.x * V2.z,
               V1.x * V2.y - V1.y * V2.x);
}

double  operator*(const Vec &V1, const Vec &V2)     // dot product
{
    return V1.x * V2.x + V1.y * V2.y + V1.z * V2.z;
}

Vec     operator+(const Vec &V, double Val)
{
    return Vec(V.x + Val, V.y + Val, V.z + Val);
}

Vec     operator-(const Vec &V, double Val)
{
    return Vec(V.x - Val, V.y - Val, V.z - Val);
}

Vec     operator*(const Vec &V, double Val)
{
    return Vec(V.x * Val, V.y * Val, V.z * Val);
}

Vec     operator/(const Vec &V, double Val)
{
    return Vec(V.x / Val, V.y / Val, V.z / Val);
}

Vec     operator+(double Val, const Vec &vector)
{
    return vector + Val;
}

Vec     operator-(double Val, const Vec &vector)
{
    return (vector - Val) * (-1.0);
}

Vec     operator*(double Val, const Vec &vector)
{
    return vector * Val;
}

Vec     operator/(double Val, const Vec &vector)
{
    return Vec(Val / vector.x,
               Val / vector.y,
               Val / vector.z);
}

// operator iostream for Vec
istream    &operator>>(istream &is, Vec &vector)
{
    is >> vector.x >> vector.y >> vector.z;
    return is;
}

ostream    &operator<<(ostream &os, const Vec &vector)
{
    os << vector.print();
    return os;
}


// ==================== Mat ====================

Mat::Mat(void)
{
    Row     = 0;
    Col     = 0;
    Data   = nullptr;
}

Mat::Mat(int Rows, int Cols, double Val)
{
    Row     = Rows;
    Col     = Cols;
    Data   = new double [size()];
    for (int i = 0; i < size(); i++)
        Data[i] = Val;
}

Mat::Mat(int Rows, int Cols, double *Array)
{
    Row     = Rows;
    Col     = Cols;
    Data   = new double [size()];
    for (int i = 0; i < size(); i++)
        Data[i] = Array[i];
}

Mat::Mat(int Rows, int Cols, double **Array)
{
    Row     = Rows;
    Col     = Cols;
    Data   = new double [size()];
    for (int i = 0; i < Row; i++)
        for (int j = 0; j < Col; j++)
            (*this)(i, j) = Array[i][j];
}

Mat::Mat(const Mat &M)
{
    Row     = M.Row;
    Col     = M.Col;
    Data   = new double [size()];
    for (int i = 0; i < size(); i++)
        Data[i] = M.Data[i];
}

Mat::~Mat(void)
{
    delete [] Data;
    Data = nullptr;
}

// operator for Mat
double &Mat::operator()(int R, int C)
{
    return Data[ R * Col + C ];
}

double  Mat::operator()(int R, int C) const
{
    return Data[ R * Col + C ];
}

Mat    &Mat::operator= (const Mat &M)
{
    if (this == &M)
        return *this;

    Row     = M.Row;
    Col     = M.Col;
    delete [] Data;
    Data   = new double [size()];
    for (int i = 0; i < size(); i++)
        Data[i] = M.Data[i];

    return *this;
}

Mat    &Mat::operator+=(const Mat &M)
{
    if (Row != M.Row || Col != M.Col)
        exit(201);

    for (int i = 0; i < size(); i++)
        Data[i] += M.Data[i];

    return *this;
}

Mat    &Mat::operator-=(const Mat &M)
{
    if (Row != M.Row || Col != M.Col)
        exit(201);

    for (int i = 0; i < size(); i++)
        Data[i] -= M.Data[i];

    return *this;
}

Mat    &Mat::operator*=(const Mat &M)
{
    if (Col != M.Row)
        exit(202);

    Mat tmp(Row, M.Col);
    for (int i = 0; i < tmp.Row; i++)
        for (int j = 0; j < tmp.Col; j++)
            for (int k = 0; k < Col; k++)
                tmp(i, j) += (*this)(i, k) * M(k, j);
    *this = tmp;

    return *this;
}

Mat    &Mat::operator+=(double Val)
{
    for (int i = 0; i < size(); i++)
        Data[i] += Val;
    return *this;
}

Mat    &Mat::operator-=(double Val)
{
    for (int i = 0; i < size(); i++)
        Data[i] -= Val;
    return *this;
}

Mat    &Mat::operator*=(double Val)
{
    for (int i = 0; i < size(); i++)
        Data[i] *= Val;
    return *this;
}

// + - * for two matrixes
Mat     operator+(const Mat &matrix)
{
    return Mat(matrix);
}

Mat     operator-(const Mat &matrix)
{
    Mat ret(matrix);
    for (int i = 0; i < ret.size(); i++)
        ret.Data[i] = - ret.Data[i];
    return ret;
}

Mat     operator+(const Mat &M1, const Mat &M2)
{
    if (M1.Row != M2.Row || M1.Col != M2.Col)
        exit(201);

    Mat ret(M1);
    for (int i = 0; i < M1.size(); i++)
        ret.Data[i] += M2.Data[i];

    return ret;
}

Mat     operator-(const Mat &M1, const Mat &M2)
{
    if (M1.Row != M2.Row || M1.Col != M2.Col)
        exit(201);

    Mat ret(M1);
    for (int i = 0; i < M1.size(); i++)
        ret.Data[i] -= M2.Data[i];

    return ret;
}

Mat     operator*(const Mat &M1, const Mat &M2)
{
    if (M1.Col != M2.Row)
        exit(202);

    Mat ret(M1.Row, M2.Col);
    for (int i = 0; i < ret.Row; i++)
        for (int j = 0; j < ret.Col; j++)
            for (int k = 0; k < M1.Col; k++)
                ret(i, j) += M1(i, k) * M2(k, j);

    return ret;
}

Mat     operator+(const Mat &M, double Val)
{
    Mat ret(M);
    for (int i = 0; i < ret.size(); i++)
        M.Data[i] += Val;

    return ret;
}

Mat     operator-(const Mat &M, double Val)
{
    Mat ret(M);
    for (int i = 0; i < ret.size(); i++)
        M.Data[i] -= Val;

    return ret;
}

Mat     operator*(const Mat &M, double Val)
{
    Mat ret(M);
    for (int i = 0; i < ret.size(); i++)
        M.Data[i] *= Val;

    return ret;
}

Mat     operator+(double Val, const Mat &M)
{
    return M + Val;
}

Mat     operator-(double Val, const Mat &M)
{
    return (M - Val) * (-1.0);
}

Mat     operator*(double Val, const Mat &M)
{
    return M * Val;
}

// operator << for Mat
ostream    &operator<<(ostream &os, const Mat &matrix)
{
    os << matrix.print();
    return os;
}

// functions of Mat

void    Mat::reshape(int Rows, int Cols)
{
    if (Rows * Cols != size())
        exit(205);

    Row = Rows;
    Col = Cols;
}

string  Mat::print() const
{
    string  ret;
    char    buff[BUFF_LEN];

    for (int i = 0; i < Row; i++) {
        for (int j = 0; j < Col; j++) {
            sprintf(buff, "%12.6lf", (*this)(i, j));
            ret += string(buff);
        }
        ret += "\n";
    }

    return ret;
}

Mat     Mat::inv()
{
    if (Row != Col)
        exit(203);

    Mat ret(*this);

    int *ipiv = new int [Row];
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, Row, Row, ret.Data, Row, ipiv);
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, Row, ret.Data, Row, ipiv);
    delete [] ipiv;
    ipiv = nullptr;

    return  ret;
}

Mat     Mat::trans()
{
    Mat ret(Col, Row);
    for (int i = 0; i < ret.Row; i++)
        for (int j = 0; j < ret.Col; j++)
            ret(i, j) = (*this)(j, i);

    return ret;
}


// ==================== EigenSolver ====================

EigenSolver::EigenSolver(const Mat &M)
{
    if (M.Row != M.Col)
        exit(203);

    Dim         = M.Row;
    EigenVecs   = Mat(Dim, Dim);
    EigenVals   = Mat(Dim, 1);

    double  *wr     = new double [Dim];
    double  *wi     = new double [Dim];
    double  *vl     = new double [M.size()];
    double  *vr     = new double [M.size()];
    double  *src    = new double [M.size()];
    for (int i = 0; i < M.size(); i++)
        src[i] = M.Data[i];

    LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', Dim, src, Dim,
                  wr, wi, vl, Dim, vr, Dim);

    // sort by eigenvalues
    int     *id = new int [Dim];
    for (int i = 0; i < Dim; i++)
        id[i] = i;
    sort(id, id + Dim, [wr](const int a, const int b) {
        return wr[a] < wr[b];
    });

    for (int i = 0; i < Dim; i++) {
        EigenVals(i, 0) = wr[id[i]];
        for (int j = 0; j < Dim; j++)
            EigenVecs(j, i) = vr[ j * Dim + id[i] ];
    }

    delete [] wr;
    wr  = nullptr;
    delete [] wi;
    wi  = nullptr;
    delete [] vl;
    vl  = nullptr;
    delete [] vr;
    vr  = nullptr;
    delete [] src;
    src = nullptr;
    delete [] id;
    id  = nullptr;
}


// calculate factor
long long   factor(int n)
{
    if (n > 20)     // overflow
        exit(101);

    long long ret = 1;
    for (int i = 1; i < n; i++)
        ret *= i;

    return ret;
}






