#include "atom.h"
#include <string>
#include <iostream>

#define BUFF_LEN 1024

using namespace std;

Atom Elem_Tab[Max_Atom_Idx + 1] = {
    Atom(0,  0.000000000, string("")),
    Atom(1,  1.007825032, string("H")),
    Atom(2,  4.002603254, string("He")),
    Atom(3,  7.016003437, string("Li")),
    Atom(4,  9.012183065, string("Be")),
    Atom(5, 11.009305364, string("B")),
    Atom(6, 12.000000000, string("C")),
    Atom(7, 14.003074004, string("N")),
    Atom(8, 15.994914620, string("O")),
    Atom(9, 18.998403163, string("F")),
    Atom(10,19.992440176, string("Ne")),
    Atom(11,0.000000000, string("Na")),
    Atom(12,0.000000000, string("Mg")),
    Atom(13,0.000000000, string("Al")),
    Atom(14,0.000000000, string("Si")),
    Atom(15,0.000000000, string("P")),
    Atom(16,0.000000000, string("S")),
    Atom(17,0.000000000, string("Cl")),
    Atom(18,0.000000000, string("Ar")),
    Atom(19,0.000000000, string("K")),
    Atom(20,0.000000000, string("Ca")),
    Atom(21,0.000000000, string("Sc")),
    Atom(22,0.000000000, string("Ti")),
    Atom(23,0.000000000, string("V")),
    Atom(24,0.000000000, string("Cr")),
    Atom(25,0.000000000, string("Mn")),
    Atom(26,0.000000000, string("Fe")),
    Atom(27,0.000000000, string("Co")),
    Atom(28,0.000000000, string("Ni")),
    Atom(29,0.000000000, string("Cu")),
    Atom(30,0.000000000, string("Zn")),
    Atom(31,0.000000000, string("Ga")),
    Atom(32,0.000000000, string("Ge")),
    Atom(33,0.000000000, string("As")),
    Atom(34,0.000000000, string("Se")),
    Atom(35,0.000000000, string("Br")),
    Atom(36,0.000000000, string("Kr"))
};


// ==================== Atom ====================

Atom::Atom(size_t idx, REAL mass, string name)
    :idx(idx), mass(mass), name(name)
{}

string  Atom::print() const
{
    char buff[BUFF_LEN];
    sprintf(buff, "%3d %-2s%14.8lf",
            INTG(idx), name.c_str(), mass);

    return  string(buff);
}

// operator << for Atom
ostream    &operator<<(ostream &os, const Atom &atom)
{
    os << atom.print();
    return os;
}


INTG    get_atom_idx(string name)
{
    INTG  ret = 0;
    for (size_t i=1; i<=Max_Atom_Idx; i++)
        if (Elem_Tab[i].name == name)
            ret = INTG(i);

    if (ret == 0) {
        cout << "no such atom:" + name << endl;
        exit(401);
    }

    return ret;
}

string  get_atom_name(INTG idx)
{
    if (size_t(idx) > Max_Atom_Idx) {
        cout << "invaild Atom Index:" + to_string(idx) << endl;
        exit(402);
    }

    return Elem_Tab[idx].name;
}
