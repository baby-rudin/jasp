#ifndef MOLECULE_H
#define MOLECULE_H

#include "mathfun.h"
#include <string>

class Molecule
{
private:
    int             nAtom;
    int             charge;
    int             *zvals;
    Vec             *geom;
    std::string     name;           // can be deleted
    std::string     point_group;

public:
    Molecule(void);
    Molecule(int nAtom);
    ~Molecule(void);

    std::string     print_geom();
    void            translate(Vec R);
    double          bond(int atom1, int atom2);
    double          angle(int atom1, int atom2, int atom3);
    double          oop(int atom1, int atom2, int atom3, int atom4);        // TO DO
    double          torsion(int atom1, int atom2, int atom3, int atom4);    // TO DO
};

#endif // MOLECULE_H
