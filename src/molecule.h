#ifndef MOLECULE_H
#define MOLECULE_H

#include "type.hpp"
#include "vector.h"
#include <string>
#include <iostream>

class Molecule
{
public:
    size_t          nAtom;
    INTG            charge;
    INTG            multipl;
    INTG            *zval;
    VecReal         *geom;
    std::string     *basis;

    Molecule(void);
    Molecule(std::istream &is);
    Molecule(const Molecule &mol);
    ~Molecule(void);

    Molecule&       operator= (const Molecule &mol);

    void            gen_basis_name(std::string basis_name);

    std::string     print_geom() const;
    void            print_info() const;
    void            translate(VecReal R);
    REAL            E_nuc() const;
    REAL            bond(size_t atom1, size_t atom2) const;
    REAL            angle(size_t atom1, size_t atom2, size_t atom3) const;
    REAL            oop( size_t atom1, size_t atom2,
                         size_t atom3, size_t atom4 ) const;        // TO DO
    REAL            torsion( size_t atom1, size_t atom2,
                             size_t atom3, size_t atom4 ) const;    // TO DO
    size_t          N_elec() const;

    // overload ostream operator
    friend std::ostream     &operator<<(std::ostream &os, const Molecule &mol);
};

#endif // MOLECULE_H
