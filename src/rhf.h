#ifndef RHF_H
#define RHF_H

#include "type.hpp"
#include "matrix.h"
#include "molecule.h"
#include "bas_set.h"

class Rhf
{
public:
    Molecule        mol;
    BasSet          bs;
    Matrix          S;          // overlap matrix
    Matrix          T;          // kinetic matrix
    Matrix          V;          // nuclear matrix
    Matrix          H;          // core Hamiltonian matrix
    Matrix          Q;          // two electron repulsion (constant)
    Matrix          X;          // transformation matrix
    Matrix          P;          // density matrix (generate by gauss first)
    Matrix          G;          // generate from P and Q
    Matrix          F;          // F = H + G
    Matrix          F_bar;      // F_bar = X' F X
    Matrix          C_bar;      // generate by diagonalize F_bar
    Matrix          E;          // generate by diagonalize F_bar
    Matrix          C;          // C = X C_bar

    REAL            E_nuc;
    REAL            E_ele;
    REAL            E_tot;


    REAL            E_old;      // used in is_converged()
    Matrix          P_old;      // used in is_converged()
    REAL            dE;         // E_tot - E_old
    REAL            dE2;

    Rhf(const Molecule &mol, const BasSet &bs);

    void    guess_P();
    void    calc_X();
    void    calc_P();
    void    calc_G();
    void    calc_F();
    void    calc_F_bar();
    void    calc_C_bar();       // calculate E simultaneously
    void    calc_C();

    void    calc_E_ele();

    bool    is_converged();

    bool    scf();
};

#endif // RHF_H
