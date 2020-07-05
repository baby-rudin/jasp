#ifndef QC_INT_H
#define QC_INT_H

#include "type.hpp"

REAL    gauss_int_overlap(REAL alpha1, INTG l1, INTG m1, INTG n1, REAL x1, REAL y1, REAL z1,
                          REAL alpha2, INTG l2, INTG m2, INTG n2, REAL x2, REAL y2, REAL z2);

REAL    gauss_int_kinetic(REAL alpha1, INTG l1, INTG m1, INTG n1, REAL x1, REAL y1, REAL z1,
                          REAL alpha2, INTG l2, INTG m2, INTG n2, REAL x2, REAL y2, REAL z2);

REAL    gauss_int_nuclear(REAL alpha1, INTG l1, INTG m1, INTG n1, REAL x1, REAL y1, REAL z1,
                          REAL alpha2, INTG l2, INTG m2, INTG n2, REAL x2, REAL y2, REAL z2,
                          REAL Zx, REAL Zy, REAL Zz);

REAL    gauss_int_repulsion(REAL alpha1, int l1, int m1, int n1, REAL x1, REAL y1, REAL z1,
                            REAL alpha2, int l2, int m2, int n2, REAL x2, REAL y2, REAL z2,
                            REAL alpha3, int l3, int m3, int n3, REAL x3, REAL y3, REAL z3,
                            REAL alpha4, int l4, int m4, int n4, REAL x4, REAL y4, REAL z4);

#endif // QC_INT_H
