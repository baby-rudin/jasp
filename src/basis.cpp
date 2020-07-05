#include "basis.h"

Basis::Basis()
    : iAtom(0), nGauss(0), gas(nullptr) {}

Basis::Basis(size_t iAtom, size_t nGauss, REAL *coeff, REAL *alpha, VecReal vec, VecIntg obt)
    : iAtom(iAtom), nGauss(nGauss), gas(nGauss > 0 ? new Gauss [nGauss] : nullptr)
{
    for (size_t i=0; i<nGauss; i++)
        gas[i] = Gauss(coeff[i], alpha[i], vec, obt);
}

Basis::Basis(const Basis &bas)
    : iAtom(bas.iAtom), nGauss(bas.nGauss), gas(bas.nGauss > 0 ? new Gauss [bas.nGauss] : nullptr)
{
    for (size_t i=0; i<nGauss; i++)
        gas[i] = bas.gas[i];
}

Basis::~Basis()
{
    delete [] gas;
    gas = nullptr;
}

Basis&  Basis::operator=(const Basis &bas)
{
    delete [] gas;

    iAtom   = bas.iAtom;
    nGauss  = bas.nGauss;
    if (nGauss > 0)
        gas = new Gauss [nGauss];
    else
        gas = nullptr;

    for (size_t i=0; i<nGauss; i++)
        gas[i] = bas.gas[i];

    return *this;
}

REAL    int_overlap(const Basis &a, const Basis &b)
{
    REAL  ret = 0.0;

    for (size_t i=0; i<a.nGauss; i++)
        for (size_t j=0; j<b.nGauss; j++)
            ret += int_overlap(a.gas[i], b.gas[j]);
    return ret;
}

REAL    int_kinetic(const Basis &a, const Basis &b)
{
    REAL  ret = 0.0;

    for (size_t i=0; i<a.nGauss; i++)
        for (size_t j=0; j<b.nGauss; j++)
            ret += int_kinetic(a.gas[i], b.gas[j]);
    return ret;
}

REAL    int_nuclear(const Basis &a, const Basis &b,     // suppose charge Zc = 1
                    const VecReal &point)
{
    REAL  ret = 0.0;

    for (size_t i=0; i<a.nGauss; i++)
        for (size_t j=0; j<b.nGauss; j++)
            ret += int_nuclear(a.gas[i], b.gas[j], point);
    return ret;
}

REAL    int_repulsion(const Basis &a, const Basis &b,
                      const Basis &c, const Basis &d)
{
    REAL  ret = 0.0;

    for (size_t i=0; i<a.nGauss; i++)
        for (size_t j=0; j<b.nGauss; j++)
            for (size_t k=0; k<c.nGauss; k++)
                for (size_t l=0; l<d.nGauss; l++)
                    ret += int_repulsion(a.gas[i], b.gas[j],
                                         c.gas[k], d.gas[l]);
    return ret;
}
