#include "gauss.h"
#include "mathfun.h"
#include "qc_int.h"
#include "constant.h"
#include "env.h"
#include <string>

#define BUFF_LEN 1024


using namespace std;

// ==================== Gauss ====================

Gauss::Gauss(REAL coeff, REAL alpha, VecReal vec, VecIntg obt)
    : coeff(coeff), vec(vec), obt(obt)
{ set_alpha(alpha); }

void    Gauss::normalize()
{
    REAL df_l = double_factorial(2 * obt.x - 1);
    REAL df_m = double_factorial(2 * obt.y - 1);
    REAL df_n = double_factorial(2 * obt.z - 1);
    norm =  pow(alpha * M_2_PI, 0.75) *
            sqrt(pow(4 * alpha, obt.sum()) /
            (df_l * df_m * df_n));
}

void    Gauss::set_alpha(REAL alpha)
{
    this->alpha = alpha;
    normalize();
}

// integral functions
REAL  int_overlap(const Gauss &a, const Gauss &b)
{
    auto vt1 = a.vec / LEN_COEFF;
    auto vt2 = b.vec / LEN_COEFF;

    return  a.norm * b.norm * a.coeff * b.coeff
          * gauss_int_overlap(a.alpha, a.obt.x, a.obt.y, a.obt.z, vt1.x, vt1.y, vt1.z,
                              b.alpha, b.obt.x, b.obt.y, b.obt.z, vt2.x, vt2.y, vt2.z);
}

REAL  int_kinetic(const Gauss &a, const Gauss &b)
{
    auto vt1 = a.vec / LEN_COEFF;
    auto vt2 = b.vec / LEN_COEFF;

    return  a.norm * b.norm * a.coeff * b.coeff
          * gauss_int_kinetic(a.alpha, a.obt.x, a.obt.y, a.obt.z, vt1.x, vt1.y, vt1.z,
                              b.alpha, b.obt.x, b.obt.y, b.obt.z, vt2.x, vt2.y, vt2.z);
}

REAL  int_nuclear(const Gauss &a, const Gauss &b, const VecReal &point)   // suppose Zc = 1 !!!
{
    auto vt1 = a.vec / LEN_COEFF;
    auto vt2 = b.vec / LEN_COEFF;
    auto mt  = point / LEN_COEFF;
    return  a.norm * b.norm * a.coeff * b.coeff
          * gauss_int_nuclear(a.alpha, a.obt.x, a.obt.y, a.obt.z, vt1.x, vt1.y, vt1.z,
                              b.alpha, b.obt.x, b.obt.y, b.obt.z, vt2.x, vt2.y, vt2.z,
                              mt.x, mt.y,mt.z);
}

REAL    int_repulsion(const Gauss &a, const Gauss &b,
                      const Gauss &c, const Gauss &d)
{
    auto vt1 = a.vec / LEN_COEFF;
    auto vt2 = b.vec / LEN_COEFF;
    auto vt3 = c.vec / LEN_COEFF;
    auto vt4 = d.vec / LEN_COEFF;
    return  a.norm * b.norm * c.norm * d.norm
          * a.coeff * b.coeff * c.coeff * d.coeff
          * gauss_int_repulsion(a.alpha, a.obt.x, a.obt.y, a.obt.z, vt1.x, vt1.y, vt1.z,
                                b.alpha, b.obt.x, b.obt.y, b.obt.z, vt2.x, vt2.y, vt2.z,
                                c.alpha, c.obt.x, c.obt.y, c.obt.z, vt3.x, vt3.y, vt3.z,
                                d.alpha, d.obt.x, d.obt.y, d.obt.z, vt4.x, vt4.y, vt4.z);
}
