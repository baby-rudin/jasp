#include "molecule.h"
#include "atom.h"
#include <string>
#include <iostream>
#include <cmath>

using namespace std;

//==================== Molecule ====================

Molecule::Molecule()
{
    nAtom   = 0;
    charge  = 0;
    zvals   = nullptr;
    geom    = nullptr;

    name        = string("");
    point_group = string("");
}

Molecule::Molecule(int nAtom)
{
    if (nAtom <= 0)
        exit(501);

    this->nAtom     = nAtom;
    this->zvals     = new int [nAtom];
    this->geom      = new Vec [nAtom];

    for (int i = 0; i < nAtom; i++) {
        string AtomSymbol;
        cin >> AtomSymbol >> geom[i];
        zvals[i] = get_atom_id(AtomSymbol);
    }
}

Molecule::~Molecule()
{
    delete [] zvals;
    zvals = nullptr;
    delete [] geom;
    geom  = nullptr;
}

// functions of Molecule
string  Molecule::print_geom()
{
    string  geometry;
    for (int i = 0; i < nAtom; i++)
        geometry +=   get_atom_symbol(zvals[i])
                      + geom[i].print()
                      + string("\n");

    return  geometry;
}

void    Molecule::translate(Vec R)
{
    for (int i = 0; i < nAtom; i++)
        geom[i] += R;
}

double  Molecule::bond(int a, int b)
{
    return (geom[a] - geom[b]).len();
}

double  Molecule::angle(int a, int b, int c)
{
    Vec b_a = geom[b] - geom[a];
    Vec b_c = geom[b] - geom[c];
    return  acos((b_a * b_c) / (b_a.len() * b_c.len()));
}



