#ifndef BASIS_H
#define BASIS_H

#include "mathfun.h"
#include <cmath>

class Gauss;    // 3d Gaussian Function
class Phi;      // 3d Math function
class Basis;

// orbit type flag
enum OrbitType {S, Px, Py, Pz};

class Gauss
{
private:
    Vec     R;
    int     ix, jx, kx;
    double  K;
    double  alpha;

    // used to compute nuclear attraction and two-electron repulsion integral
    double  F_0(double t) const
    {
        return 0.5 * pow(M_PI / t, 0.5) * erf(pow(t, 0.5));
    }

public:
    Gauss(Vec R = Vec(), double alpha = 1.0,
          int ix = 0, int jx = 0, int kx = 0);

    void    normalize();
    void    set_alpha(double alpha);

    Gauss  &operator*=(const Gauss &gaussian);

public:
    friend Gauss    operator*(const Gauss &a, const Gauss &b);

    // integral
    friend double   int_overlap(const Gauss &a, const Gauss &b);
    friend double   int_kinetic(const Gauss &a, const Gauss &b);
    friend double   int_nuclear(const Gauss &a, const Gauss &b,     // suppose charge Zc = 1
                                const Vec &point);
    friend double   int_repulsion(const Gauss &a, const Gauss &b,
                                  const Gauss &c, const Gauss &d);
};


class Phi
{
private:
    OrbitType   type;
    int         nGauss;
    Gauss       *gauss;
    double      *Di;        // coefficients
    double      zeta;

public:
    Phi(void);
    Phi(int num, double Sc, OrbitType obt,
        double *alpha, double *di, Vec point);
    ~Phi(void);

    // integral
    friend double   int_overlap(const Phi &a, const Phi &b);
    friend double   int_kinetic(const Phi &a, const Phi &b);
    friend double   int_nuclear(const Phi &a, const Phi &b,     // suppose charge Zc = 1
                                const Vec &point);
    friend double   int_repulsion(const Phi &a, const Phi &b,
                                  const Phi &c, const Phi &d);
};


class Basis
{
private:
    int     nPhi;
    Phi     *phi;

};

#endif // BASIS_H
