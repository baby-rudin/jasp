#include "basis.h"
#include <cmath>

using namespace std;

// ==================== Gauss ====================

Gauss::Gauss(Vec R, double alpha, int ix, int jx, int kx)
{
    this->R     = R;
    this->ix    = ix;
    this->jx    = jx;
    this->kx    = kx;
    set_alpha(alpha);
}

void    Gauss::normalize()
{
    double coe =  pow(8.0 * alpha, ix + jx + kx)
                  * factor(ix) * factor(jx) * factor(kx)
                  / (factor(2 * ix) * factor(2 * jx) * factor(2 * kx));

    K =  pow(2.0 * alpha / M_PI, 0.75)
         * pow(coe, 0.5);
}

void    Gauss::set_alpha(double alpha)
{
    this->alpha = alpha;
    normalize();
}

// operator *=
Gauss  &Gauss::operator*=(const Gauss &gaussian)
{
    (*this) = (*this) * gaussian;
    return *this;
}

Gauss   operator*(const Gauss &a, const Gauss &b)
{
    double  p   = a.alpha + b.alpha;
    double  Kab =  pow(2.0 * a.alpha * b.alpha / (p * M_PI), 0.75)
                   * exp(- a.alpha * b.alpha * (a.R - b.R).len2() / p);
    Vec     Rp  = (a.alpha * a.R + b.alpha * b.R) / p;
    Gauss   ret(Rp, p);
    ret.K *= Kab;

    return ret;
}

// integral functions
double  int_overlap(const Gauss &a, const Gauss &b)
{
    return a.K * b.K * pow(M_PI / (a.alpha + b.alpha), 1.5)
           * exp(-a.alpha * b.alpha / (a.alpha + b.alpha)
                 * (a.R - b.R).len2());
}

double  int_kinetic(const Gauss &a, const Gauss &b)
{
    double coeff = a.alpha * b.alpha / (a.alpha + b.alpha);
    return a.K * b.K * coeff
           * (3.0 - 2.0 * coeff * (a.R - b.R).len2())
           * pow(M_PI / (a.alpha + b.alpha), 1.5)
           * exp(-coeff * (a.R - b.R).len2());
}

double  int_nuclear(const Gauss &a, const Gauss &b, const Vec &point)   // suppose Zc = 1 !!!
{
    Gauss   p = a * b;
    return -2.0 * M_PI / (a.alpha + b.alpha)
           * exp(-a.alpha * b.alpha / (a.alpha + b.alpha) * (a.R - b.R).len2())
           * a.F_0((a.alpha + b.alpha) * (p.R - point).len2());
}

double  int_repulsion(const Gauss &a, const Gauss &b,
                      const Gauss &c, const Gauss &d)
{
    double  coeff = a.alpha + b.alpha + c.alpha + d.alpha;
    Gauss   p = a * b;
    Gauss   q = c * d;
    return a.K * b.K * c.K * d.K * 2.0 * pow(M_PI, 2.5)
           / ((a.alpha + b.alpha) * (c.alpha + d.alpha) * pow(coeff, 0.5))
           * exp(-a.alpha * b.alpha / (a.alpha + b.alpha) * (a.R - b.R).len2()
                 - c.alpha * d.alpha / (c.alpha + d.alpha) * (c.R - d.R).len2())
           * a.F_0((a.alpha + b.alpha) * (c.alpha + d.alpha) * coeff * (p.R - q.R).len2());
}


// ==================== Phi ====================

Phi::Phi(void)
{
    nGauss  = 0;
    gauss   = nullptr;
    Di      = nullptr;
    zeta    = 1.0;
}

Phi::Phi(int num, double Sc, OrbitType obt,
         double *alpha, double *di, Vec point)
{
    type    = obt;
    nGauss  = num;
    zeta    = Sc;
    gauss   = new Gauss  [nGauss];
    Di      = new double [nGauss];

    int Ix, Jx, Kx;
    switch (type) {
    case S:
        Ix = 0;
        Jx = 0;
        Kx = 0;
        break;
    case Px:
        Ix = 1;
        Jx = 0;
        Kx = 0;
        break;
    case Py:
        Ix = 0;
        Jx = 1;
        Kx = 0;
        break;
    case Pz:
        Ix = 0;
        Jx = 0;
        Kx = 1;
        break;
    }

    for (int i = 0; i < nGauss; i++) {
        Di[i] = di[i];
        gauss[i] = Gauss(point, alpha[i], Ix, Jx, Kx);
        gauss[i].set_alpha(Sc);
    }
}

Phi::~Phi(void)
{
    delete [] gauss;
    gauss = nullptr;
    delete [] Di;
    Di    = nullptr;
}


// friend function
double  int_overlap(const Phi &a, const Phi &b)
{
    double ret = 0.0;
    for (int i = 0; i < a.nGauss; i++)
        for (int j = 0; j < b.nGauss; j++)
            ret +=  a.Di[i] * b.Di[j]
                    * int_overlap(a.gauss[i], b.gauss[i]);
    return ret;
}

double  int_kinetic(const Phi &a, const Phi &b)
{
    double ret = 0.0;
    for (int i = 0; i < a.nGauss; i++)
        for (int j = 0; j < b.nGauss; j++)
            ret +=  a.Di[i] * b.Di[j]
                    * int_kinetic(a.gauss[i], b.gauss[i]);
    return ret;
}

double  int_nuclear(const Phi &a, const Phi &b,     // suppose charge Zc = 1
                    const Vec &point)
{
    double ret = 0.0;
    for (int i = 0; i < a.nGauss; i++)
        for (int j = 0; j < b.nGauss; j++)
            ret +=  a.Di[i] * b.Di[j]
                    * int_nuclear(a.gauss[i], b.gauss[i], point);
    return ret;
}

double  int_repulsion(const Phi &a, const Phi &b,
                      const Phi &c, const Phi &d)
{
    double ret = 0.0;
    for (int i = 0; i < a.nGauss; i++)
        for (int j = 0; j < b.nGauss; j++)
            for (int k = 0; k < c.nGauss; j++)
                for (int l = 0; l < d.nGauss; l++)
                    ret +=  a.Di[i] * b.Di[j] * c.Di[k] * d.Di[l]
                            * int_repulsion(a.gauss[i], b.gauss[j],
                                            c.gauss[k], d.gauss[l]);
    return ret;
}

// ==================== Basis ====================

Basis::Basis(void)
{
    nPhi    = 0;
    phi     = nullptr;
}

Basis::Basis(int)
{
    // TODO
}

Basis::~Basis(void)
{
    delete [] phi;
    phi = nullptr;
}

// friend function to calculate integral matrix;
Mat  Mat_int_overlap(const Basis &bas)      // nPhi * nPhi
{
    Mat ret(bas.nPhi, bas.nPhi);
    for (int i=0; i<bas.nPhi; i++)
        for (int j=i; j<bas.nPhi; j++)
            ret(i,j) = ret(j,i) = int_overlap(bas.phi[i], bas.phi[j]);

    return  ret;
}

Mat  Mat_int_kinetic(const Basis &bas)      // nPhi * nPhi
{
    Mat ret(bas.nPhi, bas.nPhi);
    for (int i=0; i<bas.nPhi; i++)
        for (int j=i; j<bas.nPhi; j++)
            ret(i,j) = ret(j,i) = int_kinetic(bas.phi[i], bas.phi[j]);

    return  ret;
}

Mat  Mat_int_nuclear(const Basis &bas,      // nPhi * nPhi
                     const Molecule &mol)
{
    Mat ret(bas.nPhi, bas.nPhi, 0.0);
    for (int m=0; m<mol.nAtom; m++)
        for (int i=0; i<bas.nPhi; i++)
            for (int j=i; j<bas.nPhi; j++) {
                double integral = int_nuclear(bas.phi[i], bas.phi[j], mol.geom[m])
                                  * mol.zvals[m];
                ret(i,j) += integral;
                ret(j,i) += integral;
            }

    return ret;
}

Mat  Mat_int_repulsion(const Basis &bas)       // 1 * idx(idx(nPhi), idx(nPhi))
{
    int MatSize = idx(bas.nPhi, bas.nPhi, bas.nPhi, bas.nPhi);
    Mat ret(MatSize, 1);

    for (int i=0; i<bas.nPhi; i++)
        for (int j=i; j<bas.nPhi; j++)
            for (int k=0; k<bas.nPhi; k++)
                for (int l=k; l<bas.nPhi; l++) {
                    if (idx_sym(i,j) < idx_sym(k,l))
                        continue;
                    ret(idx(i,j,k,l), 0) = int_repulsion(bas.phi[i], bas.phi[j],
                                                         bas.phi[k], bas.phi[l]);
                }

    return ret;
}

