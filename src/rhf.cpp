#include "rhf.h"
#include "env.h"
#include <iostream>

using namespace std;

Rhf::Rhf(const Molecule &mol, const BasSet &bs)
    : mol(mol), bs(bs), S(Mat_int_overlap(bs)),
      T(Mat_int_kinetic(bs)), V(Mat_int_nuclear(bs, mol)),
      H(T + V), Q(Mat_int_repulsion(bs)), E_nuc(mol.E_nuc()),
      dE(1.0)
{
    calc_X();
    guess_P();
}

void    Rhf::calc_E_ele()
{
    E_ele = 0.0;
    for (size_t i=0; i<bs.nBasis; i++)
        for (size_t j=0; j<bs.nBasis; j++)
            E_ele += P(j,i) * (F(i,j) + H(i,j));
    E_ele *= 0.5;
}


void    Rhf::guess_P()
{
    // Superposition-Of-Atomic-Densities guess
    // for the molecular density matrix

    P = P_old = Matrix(bs.nBasis, bs.nBasis);

    size_t  BsId = 0;
    INTG    zval = 0;

    for (size_t i=0; i<mol.nAtom; i++) {
        zval = mol.zval[i];

        while (zval > 0 && BsId < bs.nBasis) {
            INTG ob = bs.bas[BsId].gas[0].obt.sum();
            INTG nObt = 2 * ob + 1;

            REAL  ele_per_ob;
            if (2 * nObt < zval) {
                ele_per_ob = 2;
                zval -= 2 * nObt;
            } else {
                ele_per_ob = REAL(zval) / REAL(nObt);
                zval = 0;
            }

            for (size_t j=0; j<size_t(nObt); j++)
                P(BsId + j, BsId + j) = ele_per_ob;

            BsId += size_t(nObt);   // update basis position
        }
    }

    P /= 2.0;

}

void    Rhf::calc_P()
{
    P = Matrix(bs.nBasis, bs.nBasis);

    size_t  nElec = mol.N_elec();

    for (size_t i=0; i<bs.nBasis; i++)
        for (size_t j=0; j<bs.nBasis; j++)
            for (size_t k=0; k<nElec/2; k++)
                P(i,j) += C(i,k) * C(j,k);
    P *= 2.0;

    REAL coeff = 1.0;

    if (SCF_MIX == "AUTO") {
        if (abs(dE) < 1e-3)         coeff = 1.0;
        else if (abs(dE) < 1.0)     coeff = 0.5;
        else                        coeff = 0.1;
    }
    else if (SCF_MIX == "MIX") {
        coeff = 0.1;
    }
    else {
        coeff = 1.0;
    }

    P = coeff * P + (1.0 - coeff) * P_old;
}

void    Rhf::calc_G()
{
    G = Matrix(bs.nBasis, bs.nBasis);

    for (size_t i=0; i<bs.nBasis; i++)
    for (size_t j=0; j<bs.nBasis; j++)
        for (size_t k=0; k<bs.nBasis; k++)
        for (size_t l=0; l<bs.nBasis; l++)
                    G(i,j) += P(k,l)
                            *(
                               Q(idx4(i,j,l,k)) -
                               0.5 * Q(idx4(i,k,l,j))
                             );
}

void    Rhf::calc_F()
{
    F = H + G;
}

void    Rhf::calc_X()
{
    EigenSolver solver(S);

    Matrix U   = solver.eig_vec();
    Matrix s   = solver.eig_val();

    Matrix s_1_2(s);

    for (size_t i=0; i<bs.nBasis; i++)
        s_1_2(i,i) = 1.0 / sqrt(s_1_2(i,i));

    X = U * s_1_2 * U.trans();
}

void    Rhf::calc_F_bar()
{
    F_bar = X.trans() * F * X;
}

void    Rhf::calc_C_bar()
{
    EigenSolver solver(F_bar);
    C_bar   = solver.eig_vec();
    E       = solver.eig_val();
}

void    Rhf::calc_C()
{
    C = X * C_bar;
}

bool    Rhf::is_converged()
{
    bool    ret = false;
    size_t  K   = bs.nBasis;
    Matrix  tmp = P - P_old;
    REAL    err = 0.0;


    for (size_t i=0; i<bs.nBasis; i++)
        for (size_t j=0; j<bs.nBasis; j++)
            err += tmp(i,j) * tmp(i,j);
    dE2 = dE;
    dE  = E_tot - E_old;
    if ( sqrt(err/(K*K)) < SCF_P_EPS &&
         abs(dE) < SCF_E_EPS ) {
        ret = true;
    }

    // update P_old
    P_old = P;
    E_old = E_tot;

    return ret;
}


bool    Rhf::scf()
{
    if (PRINT_LEVEL == "ALL" ||
        PRINT_LEVEL == "MORE" ) {
        cout << "*** Overlap *** " << endl;
        cout << S << endl;
        cout << "*** Kinetic Energy ***" << endl;
        cout << T << endl;
        cout << "***** Potential Energy *****" << endl;
        cout << V << endl;
        cout << "****** Core Hamiltonian ******" << endl;
        cout << H << endl;
    }

    if (PRINT_LEVEL == "ALL") {
        cout << "****** Repulsion Integral ******" << endl;
        size_t  nBas = bs.nBasis;
        for (size_t i=0; i<nBas; i++)
        for (size_t j=i; j<nBas; j++) {
            size_t ij = idx2(i,j);
            for (size_t k=0; k<nBas; k++)
            for (size_t l=k; l<nBas; l++) {
                size_t kl = idx2(k,l);
                if (ij >= kl)
                    printf(" i =%3d     j =%3d     k =%3d     "
                           "l = %3d      value =  %9.6lf\n",
                           INTG(i+1), INTG(j+1), INTG(k+1), INTG(l+1),
                           Q(idx2(ij,kl)));
            }
        }
        cout << endl;
    }

    bool    is_conv = false;
    INTG    cycles  = 0;

    while (true) {
        if ( (++cycles) > SCF_MAX_ROUND ) {
            cout << "scf is Not converged after " << SCF_MAX_ROUND << " cycles  " << endl;
            break;
        }

        calc_G();
        calc_F();
        calc_F_bar();
        calc_C_bar();
        calc_C();
        calc_P();

        calc_E_ele();
        E_tot = E_ele + E_nuc;

        is_conv = is_converged();

        if (PRINT_LEVEL != "LESS" ) {

            cout << "                  cycle:       " << cycles << endl;
            cout << " ************************************************** " << endl;

            printf("  E_tot = %16.10lf", E_tot);
            if (cycles == 1) printf("\n");
            else printf("      dE = %12.10lf\n", dE);
            printf("  E_nuc = %16.10lf\n", E_nuc);
            printf("  E_ele = %16.10lf\n", E_ele);

            cout << "\n" << endl;
        }

        if (is_conv) break;
    }

    if (is_conv) {
        printf(" SCF Done: E(RHF) = %16.10lf     A.U. after %4d cycles\n",
               E_tot, cycles);
        printf("           E_ele = %16.10lf     E_nuc = %16.10lf\n",
               E_ele, E_nuc);

        cout << endl;
    }

    return is_conv;
}
