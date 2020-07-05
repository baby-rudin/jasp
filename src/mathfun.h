#ifndef MATHFUN_H
#define MATHFUN_H

#include "type.hpp"
#include "constant.h"

REAL    factorial(INTG n);          // factorial: n!
REAL    double_factorial(INTG n);   // double factorial: n!!
REAL    combntns(INTG n, INTG x);   // Combination number:  C_n^x
REAL    F_m(INTG m, REAL t);        // incomplete gamma function

// the coeffcient of x^j  in the expansion of (x+a)^l * (x+b)^m
REAL    binom_coeff(INTG j, INTG l, INTG m, REAL a, REAL b);

#endif // MATHFUN_H
