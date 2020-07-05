#ifndef GAUSS_H
#define GAUSS_H

#include "vector.h"
#include "type.hpp"
#include <string>

class Gauss
{
public:
    REAL      coeff;      // coefficient for combinating

    REAL      alpha;
    REAL      norm;
    VecReal   vec;
    VecIntg   obt;

    Gauss(REAL coeff = 1.0, REAL alpha = 1.0,
          VecReal vec = VecReal(), VecIntg obt = VecIntg());

    void            normalize();
    void            set_alpha(REAL alpha);

    // integral
    friend REAL   int_overlap(const Gauss &a, const Gauss &b);
    friend REAL   int_kinetic(const Gauss &a, const Gauss &b);
    friend REAL   int_nuclear(const Gauss &a, const Gauss &b,     // suppose charge Zc = 1
                              const VecReal &point);
    friend REAL   int_repulsion(const Gauss &a, const Gauss &b,
                                const Gauss &c, const Gauss &d);

};

REAL   int_overlap(const Gauss &a, const Gauss &b);
REAL   int_kinetic(const Gauss &a, const Gauss &b);
REAL   int_nuclear(const Gauss &a, const Gauss &b,     // suppose charge Zc = 1
                   const VecReal &point);
REAL   int_repulsion(const Gauss &a, const Gauss &b,
                     const Gauss &c, const Gauss &d);

#endif // GAUSS_H
