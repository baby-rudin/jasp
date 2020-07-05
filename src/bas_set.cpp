#include "bas_set.h"
#include "atom.h"
#include "stropt.h"
#include "basis.h"
#include "mathfun.h"
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

BasSet::BasSet()
    : nBasis(0), bas(nullptr)
{}

BasSet::BasSet(const Molecule &mol)
{
    vector<Basis> BasTmp;

    for (size_t i=0; i<mol.nAtom; i++) {
        ifstream is("basis/" + mol.basis[i]);

        if (!is.good()) {
            cout << "Can't open file: basis/" + mol.basis[i] << endl;
            exit(1353);
        }

        string line;
        while (getline(is, line)) {
            clean_line(line);
            if (line.substr(0,4) == "****")
                break;
        }

        while (getline(is, line)) {
            clean_line(line);
            stringstream ss(line);
            string  symbol;
            int     zero;
            ss >> symbol >> zero;

            while (getline(is, line)) {
                clean_line(line);

                if (line.substr(0,4) == "****")
                    break;

                if (symbol == get_atom_name(mol.zval[i])) {
                    string  name;
                    size_t  nGauss;
                    REAL  Sc;
                    stringstream ss1(line);
                    ss1 >> name >> nGauss >> Sc;
                    REAL  *alpha = new REAL [nGauss];
                    REAL  *comb1 = new REAL [nGauss];
                    REAL  *comb2 = new REAL [nGauss];

                    for (size_t j=0; j<nGauss; j++) {
                        string  line2;
                        getline(is,line2);
                        clean_line(line2);
                        str_change(line, 'D', 'E');     // for c++ stream
                        stringstream ss2(line2);
                        ss2 >> alpha[j] >> comb1[j];
                        if (name == "SP")
                            ss2 >> comb2[j];
                    }

                    if (name == "S") {
                        BasTmp.push_back(Basis(i, nGauss, comb1, alpha, mol.geom[i], VecIntg(0,0,0)));
                    }
                    else if (name == "SP") {
                        BasTmp.push_back(Basis(i, nGauss, comb1, alpha, mol.geom[i], VecIntg(0,0,0)));
                        BasTmp.push_back(Basis(i, nGauss, comb2, alpha, mol.geom[i], VecIntg(1,0,0)));
                        BasTmp.push_back(Basis(i, nGauss, comb2, alpha, mol.geom[i], VecIntg(0,1,0)));
                        BasTmp.push_back(Basis(i, nGauss, comb2, alpha, mol.geom[i], VecIntg(0,0,1)));
                    }
                    else if (name == "D") {
                        BasTmp.push_back(Basis(i, nGauss, comb1, alpha, mol.geom[i], VecIntg(2,0,0)));
                        BasTmp.push_back(Basis(i, nGauss, comb1, alpha, mol.geom[i], VecIntg(0,2,0)));
                        BasTmp.push_back(Basis(i, nGauss, comb1, alpha, mol.geom[i], VecIntg(0,0,2)));
                        BasTmp.push_back(Basis(i, nGauss, comb1, alpha, mol.geom[i], VecIntg(0,1,1)));
                        BasTmp.push_back(Basis(i, nGauss, comb1, alpha, mol.geom[i], VecIntg(1,0,1)));
                        BasTmp.push_back(Basis(i, nGauss, comb1, alpha, mol.geom[i], VecIntg(1,1,0)));
                    }

                    delete [] alpha;    alpha = nullptr;
                    delete [] comb1;    comb1 = nullptr;
                    delete [] comb2;    comb2 = nullptr;
                }
            }
        }
        is.close();
    }

    nBasis = BasTmp.size();
    bas    = new Basis [nBasis];
    for (size_t i=0; i<nBasis; i++)
        bas[i] = BasTmp[i];
}

BasSet::BasSet(const BasSet &bSet)
{
    nBasis  = bSet.nBasis;
    bas     = new Basis [nBasis];
    for (size_t i=0; i<nBasis; i++)
        bas[i] = bSet.bas[i];
}

BasSet::~BasSet()
{
    delete [] bas;
    bas = nullptr;
}

BasSet  &BasSet::operator=(const BasSet &bSet)
{
    delete [] bas;

    nBasis  = bSet.nBasis;
    bas     = new Basis [nBasis];
    for (size_t i=0; i<nBasis; i++)
        bas[i] = bSet.bas[i];

    return *this;
}

Matrix  Mat_int_overlap(const BasSet &bs)
{
    cout << "Calculating Overlap Integrals ..." << endl;

    Matrix  ret(bs.nBasis, bs.nBasis);

    for (size_t i=0; i<bs.nBasis; i++)
        for (size_t j=i; j<bs.nBasis; j++)
            ret(i,j) = ret(j,i) = int_overlap(bs.bas[i], bs.bas[j]);
    cout << "Done.\n" << endl;

    return ret;
}

Matrix   Mat_int_kinetic(const BasSet &bs)
{
    cout << "Calculating Kinetic Integrals ..." << endl;

    Matrix  ret(bs.nBasis, bs.nBasis);

    for (size_t i=0; i<bs.nBasis; i++)
        for (size_t j=i; j<bs.nBasis; j++)
            ret(i,j) = ret(j,i) = int_kinetic(bs.bas[i], bs.bas[j]);
    cout << "Done.\n" << endl;

    return ret;
}

Matrix   Mat_int_nuclear(const BasSet &bSet,     // nPhi * nPhi
                         const Molecule &mol)
{
    cout << "Calculating Nuclear Integrals ..." << endl;

    Matrix  ret(bSet.nBasis, bSet.nBasis);

    for (size_t i=0; i<bSet.nBasis; i++) {
        for (size_t j=i; j<bSet.nBasis; j++) {
            REAL sum = 0.0;
            for (size_t k=0; k<mol.nAtom; k++)
                sum += int_nuclear(bSet.bas[i], bSet.bas[j], mol.geom[k]) * mol.zval[k];
            ret(i,j) = ret(j,i) = sum;
        }
    }
    cout << "Done.\n" << endl;

    return ret;
}

Matrix   Mat_int_repulsion(const BasSet &bs)  // 1 * idx(idx(nPhi), idx(nPhi))
{
    cout << "Calculating Repulsion Integrals ..." << endl;

    size_t  size = idx4(bs.nBasis-1, bs.nBasis-1,
                        bs.nBasis-1, bs.nBasis-1) + 1;

    Matrix  ret(size);

    for (size_t i=0; i<bs.nBasis; i++) {
    for (size_t j=0; j<=i; j++) {
        size_t ij = idx2(i,j);
        for (size_t k=0; k<bs.nBasis; k++) {
        for (size_t l=0; l<=k; l++) {
            size_t  kl = idx2(k,l);
            if (ij <= kl)
                ret(idx2(ij, kl)) = int_repulsion(bs.bas[i], bs.bas[j],
                                                  bs.bas[k], bs.bas[l]);
        }
        }
    }
    }
    cout << "Done.\n" << endl;

    return  ret;
}


// used for calculating position of repulsion integral
// in a one-dim matrix
size_t  idx2(size_t i, size_t j)
{ return  i>j ? i * (i+1) / 2 + j : j * (j+1) / 2 +i; }

size_t  idx4(size_t i, size_t j, size_t k, size_t l)
{ return  idx2(idx2(i,j), idx2(k,l)); }
