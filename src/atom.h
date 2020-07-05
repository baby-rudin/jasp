#ifndef ATOM_H
#define ATOM_H

#include "type.hpp"
#include <string>

class Atom
{
public:
    size_t          idx;
    REAL            mass;
    std::string     name;

    Atom(size_t idx = 0, REAL mass = 0.0,
         std::string name = std::string());

    std::string     print() const;

    friend std::ostream &operator<<(std::ostream &os, const Atom &atom);
};

const size_t    Max_Atom_Idx = 36;
extern Atom     Elem_Tab[Max_Atom_Idx + 1];

INTG            get_atom_idx(std::string name);
std::string     get_atom_name(INTG idx);

#endif // ATOM_H
