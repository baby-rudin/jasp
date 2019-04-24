#ifndef SCF_H
#define SCF_H

#include "mathfun.h"
#include "basis.h"
#include <string>

extern double   scf_P_eps;
extern int      scf_max_round;

class Roothaan
{
private:
    Molecule    mol;
    Basis       bas;
    Mat         S;          // overlap matrix
    Mat         H;          // core Hamiltonian matrix
    Mat         Q;          // two electron repulsion (constant)
    Mat         X;          // transformation matrix
    Mat         P;          // density matrix (generate by gauss first)
    Mat         G;          // generate from P and Q
    Mat         F;          // F = H + G
    Mat         F_bar;      // F_bar = X' F X
    Mat         C_bar;      // generate by diagonalize F_bar
    Mat         E;          // generate by diagonalize F_bar
    Mat         C;          // C = X C_bar
    Mat         P_old;      // used in is_converged()

public:
    Roothaan(const Molecule &molecule, const Basis &basis);

    void    guess_P();
    void    calc_P(int nElec);
    void    calc_G();
    void    calc_F();
    void    calc_F_bar();
    void    calc_C_bar();       // calculate E simultaneously
    void    calc_C();
    bool    is_converged();

    std::string     print_info(int level);

    friend void     scf(Roothaan &roothaan, int print);
};

#endif // SCF_H
