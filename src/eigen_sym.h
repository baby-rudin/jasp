#ifndef EIGEN_SYM_H
#define EIGEN_SYM_H

#include "type.hpp"

// solve eigen problem of symmetric matrix
bool eig_sym(REAL *pMatrix, INTG nDim, REAL *eigVal, REAL *eigVec, REAL dbEps = 1e-10,  INTG nJt = -1);

#endif // EIGEN_SYM_H
