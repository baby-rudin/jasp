#include "atom.h"
#include <string>
#include <iostream>

#define BUFF_LEN 1024

using namespace std;

Atom Elem_Tab[] = {
    Atom(0,  0.000000000, string("")),
    Atom(1,  1.007825032, string("H")),
    Atom(2,  4.002603254, string("He")),
    Atom(3,  7.016003437, string("Li")),
    Atom(4,  9.012183065, string("Be")),
    Atom(5, 11.009305364, string("B"))
};


// ==================== Atom ====================

Atom::Atom(void)
{
    id      = 0;
    mass    = 0.0;
    symbol  = string("");
}

Atom::Atom(int id, double mass, string symbol)
{
    this->id        = id;
    this->mass      = mass;
    this->symbol    = symbol;
}

string  Atom::print() const
{
    char buff[BUFF_LEN];
    sprintf(buff, "%3d  %-2s%14.8lf",
            id, symbol.c_str(), mass);

    return  string(buff);
}

// operator << for Atom
ostream    &operator<<(ostream &os, const Atom &atom)
{
    os << atom.print();
    return os;
}


int     get_atom_id(string symbol)
{
    int ret = -1;
    for (int i = 1; i <= Max_Atom_Id; i++)
        if (Elem_Tab[i].symb() == symbol)
            ret = i;

    if (ret == -1)
        exit(401);

    return ret;
}

string  get_atom_symbol(int id)
{
    if (id <= 0 || id > Max_Atom_Id)
        exit(402);

    return Elem_Tab[id].symb();
}
