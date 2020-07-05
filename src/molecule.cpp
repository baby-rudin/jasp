#include "molecule.h"
#include "constant.h"
#include "atom.h"
#include "stropt.h"
#include "env.h"
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#define BUFF_LEN 1024

using namespace std;

//==================== Molecule ====================
Molecule::Molecule(void) {}

Molecule::Molecule(istream  &is)
{
    string  line;

    // read in charge and multi
    while (getline(is, line)) {
        clean_line(line);
        if (!line.empty())
            break;
    }


    stringstream    ss(line);
    ss >> charge >> multipl;

    // read in geometry
    vector<INTG>        ZvalTmp;
    vector<VecReal>     GeomTmp;
    while (getline(is, line)) {
        clean_line(line);
        if (line.empty())
            break;

        stringstream    ssl(line);
        string  symbol;
        VecReal position;
        ssl >> symbol >> position.x >> position.y >> position.z;
        ZvalTmp.push_back(INTG(get_atom_idx(symbol)));
        GeomTmp.push_back(position);
    }

    nAtom   = ZvalTmp.size();
    zval    = new INTG    [nAtom];
    geom    = new VecReal [nAtom];
    basis   = new string  [nAtom];
    for (size_t i=0; i<nAtom; i++) {
        zval[i] = ZvalTmp[i];
        geom[i] = GeomTmp[i];
    }
}

Molecule::Molecule(const Molecule &mol)
    : nAtom(mol.nAtom), charge(mol.charge), multipl(mol.multipl),
      zval(nAtom ? new INTG [nAtom] : nullptr), 
      geom(nAtom ? new VecReal [nAtom] : nullptr),
      basis(nAtom ? new string [nAtom] : nullptr)
{
    for (size_t i=0; i<nAtom; i++) {
        zval[i] = mol.zval[i];
        geom[i] = mol.geom[i];
        basis[i]= mol.basis[i];
    }
}

Molecule::~Molecule()
{
    delete [] zval;     zval  = nullptr;
    delete [] geom;     geom  = nullptr;
    delete [] basis;    basis = nullptr;
}

void    Molecule::gen_basis_name(std::string basis_name)
{
    str_upper(basis_name);
    if (basis_name.substr(0,3) == "GEN") {
        string line;
        while (getline(cin, line)) {
            clean_line(line);
            if (!line.empty())
                break;
        }

        while (true) {
            if (!line.empty())
                break;

            string  _basis;
            getline(cin, _basis);

            stringstream ss(line);
            string symb;
            while( ss >> symb && symb != "0" ) {
                INTG zid = get_atom_idx(symb);
                for (size_t i=0; i<nAtom; i++)
                    if (zval[i] == zid)
                        basis[i] = _basis;
            }

            getline(cin, line);
            getline(cin, line);
            clean_line(line);
        }
    }
    else {
        for (size_t i=0; i<nAtom; i++)
            basis[i] = basis_name;
    }
}

// functions of Molecule
string  Molecule::print_geom() const
{
    string  geometry;
    char buff[BUFF_LEN];

    for (size_t i = 0; i < nAtom; i++) {
        sprintf(buff, " %-2s %14.8lf%14.8lf%14.8lf",
                get_atom_name(zval[i]).c_str(), geom[i].x, geom[i].y, geom[i].z);
        geometry += string(buff) + string("\n");
    }
    return  geometry;
}

void    Molecule::print_info() const
{
    cout << " **** Charge and Multipl ****" << endl;
    cout << " Charge = " << charge << "      "
         << "Multipl = " << multipl << endl;
    cout << endl;

    string  unit_str;
    if (UNIT_LENGTH == "ANG") unit_str = "Angstrom";
    else if (UNIT_LENGTH == "BOHR") unit_str = "Bohr";

    cout << "                    Molecule orientation:                      " << endl;
    cout << " ------------------------------------------------------------- " << endl;
    cout << " Center     Atomic               Coordinates (" << unit_str << ")       " << endl;
    cout << " Number     Number           X             Y             Z     " << endl;
    cout << " ------------------------------------------------------------- " << endl;
    for (size_t i=0; i<nAtom; i++) {
        printf(" %3d         %2d     ", INTG(i+1), zval[i]);
        printf("%14.8lf%14.8lf%14.8lf\n", geom[i].x, geom[i].y, geom[i].z);
    }
    cout << " ------------------------------------------------------------- \n" << endl;
}

void    Molecule::translate(VecReal R)
{
    for (size_t i = 0; i < nAtom; i++)
        geom[i] += R;
}

REAL  Molecule::E_nuc() const
{
    REAL  ret = 0.0;
    for (size_t i=0; i<nAtom; i++)
        for (size_t j=i+1; j<nAtom; j++)
            ret += zval[i] * zval[j] * ANGSTROM_PER_BOHR / bond(i,j);
    return ret;
}

REAL  Molecule::bond(size_t a, size_t b) const
{
    return (geom[a] - geom[b]).len();
}

REAL  Molecule::angle(size_t a, size_t b, size_t c) const
{
    VecReal b_a = geom[b] - geom[a];
    VecReal b_c = geom[b] - geom[c];
    return  acos((b_a * b_c) / (b_a.len() * b_c.len()));
}

size_t  Molecule::N_elec() const
{
    int sum = 0;
    for (size_t i = 0; i<nAtom; i++)
        sum += zval[i];
    sum -= charge;
    return size_t(sum);
}

// operator =
Molecule&   Molecule::operator= (const Molecule &mol)
{
    delete [] zval;     zval  = nullptr;
    delete [] geom;     geom  = nullptr;
    delete [] basis;    basis = nullptr;

    charge  = mol.charge;
    multipl = mol.multipl;
    nAtom   = mol.nAtom;
    zval    = new INTG    [nAtom];
    geom    = new VecReal [nAtom];
    basis   = new string  [nAtom];
    for (size_t i=0; i<nAtom; i++) {
        zval[i] = mol.zval[i];
        geom[i] = mol.geom[i];
        basis[i]= mol.basis[i];
    }

    return *this;
}

ostream&    operator<<(std::ostream &os, const Molecule &mol)
{
    os << mol.charge << " " << mol.multipl << endl;
    os << mol.print_geom();
    return os;
}
