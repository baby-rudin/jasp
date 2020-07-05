#include "mathfun.h"
#include <cmath>
#include <iostream>

using namespace std;

REAL    factorial(INTG n)
{
    if (n < 0) {
        cout << "negative factorial!" << endl;
        exit(2324);
    }


    REAL ret = 1.0;

    switch (n) {
        case 0: ret = 1.0; break;
        case 1: ret = 1.0; break;
        case 2: ret = 2.0; break;
        case 3: ret = 6.0; break;
        case 4: ret = 24.0; break;
        case 5: ret = 120.0; break;
        case 6: ret = 720.0; break;
        case 7: ret = 5040.0; break;
        case 8: ret = 40320.0; break;
        case 9: ret = 362880.0; break;
        case 10: ret = 3628800.0; break;
        case 11: ret = 39916800.0; break;
        case 12: ret = 479001600.0; break;
        case 13: ret = 6227020800.0; break;
        case 14: ret = 87178291200.0; break;
        case 15: ret = 1307674368000.0; break;
        case 16: ret = 20922789888000.0; break;
        case 17: ret = 355687428096000.0; break;
        case 18: ret = 6402373705728000.0; break;
        case 19: ret = 121645100408832000.0; break;
        case 20: ret = 2432902008176640000.0; break;
        default:
            for (INTG i=1; i<=n; i++)
                ret *= i;
    }

    return ret;
}

REAL    double_factorial(INTG n)
{
    REAL ret = 1.0;

    if (n >= 0) {
        for (INTG i=((n&1) ? 1 : 2); i<=n; i+=2)
            ret *= i;
    }
    else if (n&1) {
        for (INTG i=1; i<abs(n); i+=2)
            ret /= i;
    }
    else {
        cout << "negative even double factorial!" << endl;
        exit(15616);
    }
    return  ret;
}

REAL    combntns(INTG n, INTG x)
{
    REAL ret = 1.0;

    if (x > n-x)  x = n - x;

    for(INTG i=n; i>n-x; i--)
        ret *= i;
    for(INTG i=1; i<=x; i++)
        ret /= i;

    return ret;
}

REAL    F_m(INTG m, REAL t)
{
    if (t < 1e-5)
        return 1.0 / (2.0 * m + 1);
    else if (m == 0)
        return 0.5 * sqrt(M_PI/t) * erf(sqrt(t));
    else
        return 0.5 / t * ( (2*m -1) * F_m(m - 1, t) - exp(-t) );
}

REAL    binom_coeff(INTG j, INTG l, INTG m, REAL a, REAL b)
{
    REAL ret = 0.0;

    for (INTG p=0; p<=l; p++) {
        INTG q = j-p;
        if (q >= 0 && q <= m)
            ret += combntns(l, p) *combntns(m, q) *
                   pow(a, l-p) * pow(b, m-q);
    }

    return ret;
}

