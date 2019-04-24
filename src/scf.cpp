#include "scf.h"
#include "mathfun.h"
#include "basis.h"
#include "molecule.h"
#include "atom.h"
#include <string>
#include <cstdio>
#include <iostream>

using namespace std;

double  scf_P_eps       = 1e-4;
int     scf_max_round   = 300;

// ==================== Roothaan ====================

Roothaan::Roothaan(const Molecule &molecule, const Basis &basis)
{
    mol     = molecule;
    bas     = basis;
    S       = Mat_int_overlap(bas);
    H       = Mat_int_kinetic(bas) + Mat_int_nuclear(bas, mol);
    Q       = Mat_int_repulsion(bas);
    guess_P();
}

void    Roothaan::guess_P()
{
    P = P_old = Mat(bas.size(), bas.size());
}

void    Roothaan::calc_P(int nElec)
{
    P = Mat(bas.size(), bas.size());

    for (int i=0; i<bas.size(); i++)
        for (int j=0; j<bas.size(); j++)
            for (int k=0; k<nElec/2; k++)
                P(i,j) += C(i,k) * C(j,k);
    P *= 2.0;
}

void    Roothaan::calc_G()
{
    G = Mat(bas.size(), bas.size());

    for (int i=0; i<bas.size(); i++)
        for (int j=0; j<bas.size(); j++)
            for (int k=0; k<bas.size(); k++)
                for (int l=0; l<bas.size(); l++)
                    G(i,j) += P(k,l)
                              *(Q(idx(i,j,l,k),0) - 0.5 * Q(idx(i,k,l,j),0));
}

void    Roothaan::calc_F()
{
    F = H + G;
}

void    Roothaan::calc_F_bar()
{
    EigenSolver solver(S);

    Mat U   = solver.eig_vec();
    Mat s   = solver.eig_val();

    Mat s_1_2(s);

    for (int i=0; i<bas.size(); i++)
        s_1_2(i,i) = sqrt(s_1_2(i,i));

    X = U * s_1_2 * U.trans();

    F_bar = X.trans() * F * X;
}

void    Roothaan::calc_C_bar()
{
    EigenSolver solver(F_bar);
    C_bar   = solver.eig_vec();
    E       = solver.eig_val();
}

void    Roothaan::calc_C()
{
    C = X * C_bar;
}

bool    Roothaan::is_converged()
{
    bool    ret = false;
    Mat     tmp = P - P_old;
    double  sum = 0.0;

    for (int i=0; i<bas.size(); i++)
        for (int j=0; j<bas.size(); j++)
            sum += tmp(i,j) * tmp(i,j);

    if (pow(sum/(bas.size()*bas.size()), 0.5) < scf_P_eps)
        ret = true;

    // update P_old
    P_old = P;

    return ret;
}

string  Roothaan::print_info(int level)
{
    if (level) {
        cout << "C:" << endl;
        cout << C << endl;
        cout << "P:" << endl;
        cout << P << endl;
        cout << "G:" << endl;
        cout << G << endl;
        cout << "F:" << endl;
        cout << F << endl;
        cout << "F_bar:" << endl;
        cout << F_bar << endl;
        cout << "C_bar:" << endl;
        cout << C_bar << endl;
        cout << "E:" << endl;
        cout << E << endl;
    }
}



// ========== scf function ==========
void    scf(Roothaan &rooth, int print)
{
    if (print) {
        cout << "H:" << endl;
        cout << rooth.H << endl;
        cout << "S:" << endl;
        cout << rooth.S << endl;
        cout << "X:" << endl;
        cout << rooth.X << endl;
    }

    int round = 0;

    do{
        round ++;
        rooth.calc_P();
        rooth.calc_G();
        rooth.calc_F();
        rooth.calc_F_bar();
        rooth.calc_C_bar();
        rooth.calc_C();
        rooth.print_info(print);
    } while (!rooth.is_converged() && round < scf_max_round);
}












