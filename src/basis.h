#ifndef BASIS_H
#define BASIS_H

#include "gauss.h"
#include "vector.h"
#include "type.hpp"


class Basis
{
public:
    size_t      iAtom;      // which atom it belongs to
    size_t      nGauss;
    Gauss       *gas;

    Basis();
    Basis(size_t iAtom, size_t nGauss, REAL *coeff, REAL *alpha, VecReal vec, VecIntg obt);
    Basis(const Basis &bas);
    ~Basis();

    Basis   &operator=(const Basis &bas);

    friend REAL     int_overlap(const Basis &a, const Basis &b);
    friend REAL     int_kinetic(const Basis &a, const Basis &b);
    friend REAL     int_nuclear(const Basis &a, const Basis &b,     // suppose charge Zc = 1
                                const VecReal &point);
    friend REAL     int_repulsion(const Basis &a, const Basis &b,
                                  const Basis &c, const Basis &d);
};

REAL    int_overlap(const Basis &a, const Basis &b);
REAL    int_kinetic(const Basis &a, const Basis &b);
REAL    int_nuclear(const Basis &a, const Basis &b,     // suppose charge Zc = 1
                    const VecReal &point);
REAL    int_repulsion(const Basis &a, const Basis &b,
                      const Basis &c, const Basis &d);

#endif // BASIS_H
