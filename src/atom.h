#ifndef ATOM_H
#define ATOM_H

#include <string>

const int       Max_Atom_Id = 5;

class Atom;
extern Atom     Elem_Tab[Max_Atom_Id + 1];

class Atom
{
private:
    int             id;
    double          mass;
    std::string     symbol;

public:
    Atom(void);
    Atom(int id, double mass, std::string symbol);

    std::string     print() const;
    std::string     symb() const
    {
        return symbol;
    }

public:
    friend std::ostream &operator<<(std::ostream &os, const Atom &atom);
};


int             get_atom_id(std::string symbol);
std::string     get_atom_symbol(int id);


#endif // ATOM_H
