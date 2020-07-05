#ifndef BAS_SET_H
#define BAS_SET_H


#include "basis.h"
#include "molecule.h"
#include "matrix.h"
#include <string>

class BasSet
{
public:
    size_t      nBasis;
    Basis       *bas;

    BasSet();
    BasSet(const Molecule &mol);
    BasSet(const BasSet &bSet);
    ~BasSet();

    BasSet  &operator=(const BasSet &bSet);

    friend Matrix   Mat_int_overlap(const BasSet &bs);    // nPhi * nPhi
    friend Matrix   Mat_int_kinetic(const BasSet &bs);    // nPhi * nPhi
    friend Matrix   Mat_int_nuclear(const BasSet &bs,     // nPhi * nPhi
                                    const Molecule &mol);
    friend Matrix   Mat_int_repulsion(const BasSet &bs);  // 1 * idx(idx(nPhi), idx(nPhi))
};

Matrix  Mat_int_overlap(const BasSet &bs);    // nPhi * nPhi
Matrix  Mat_int_kinetic(const BasSet &bs);    // nPhi * nPhi
Matrix  Mat_int_nuclear(const BasSet &bs,     // nPhi * nPhi
                        const Molecule &mol);
Matrix  Mat_int_repulsion(const BasSet &bs);  // 1 * idx(idx(nPhi), idx(nPhi))


// used for calculating position of repulsion integral
// in a one-dim matrix
size_t  idx2(size_t i, size_t j);
size_t  idx4(size_t i, size_t j, size_t k, size_t l);

#endif // BAS_SET_H
