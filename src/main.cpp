#include <iostream>
#include <chrono>
#include <ctime>
#include "env.h"
#include "vector.h"
#include "molecule.h"
#include "bas_set.h"
#include "rhf.h"

#define BUFF_LEN 1024

using namespace std;

int main()
{
    char buff[BUFF_LEN];

    // output head
    auto start = std::chrono::system_clock::now();
    auto start_time = std::chrono::system_clock::to_time_t(start);

    cout << " ****************************************** " << endl;
    cout << "      Jasp:  " << ctime(&start_time);
    cout << " ****************************************** " << endl;
    cout << " " << endl;

    // read in commands and
    // set Environmental Variables
    string  command_line;
    getline(cin, command_line);
    if (!set_env(command_line)) {
        cout << "Error control line!" << endl;
        exit(127);
    }

    // output control info
    print_ctrl_info();

    // read in molecule
    Molecule mol(cin);
    mol.gen_basis_name(BASIS_SETS);

    // print input molecule information
    mol.print_info();

    // read in basis set
    BasSet  bs(mol);

    // calculate
    bool is_normal = false;
    if ( CALC_METHOD == "RHF" ||
         CALC_METHOD == "HF"     ) {
        Rhf rhf(mol, bs);
        is_normal = rhf.scf();
    }
    else if (CALC_METHOD == "UHF") {
        cout << "Method " + CALC_METHOD + "haven't supported now!" << endl;
    }
    else {
        cout << "Method " + CALC_METHOD + "haven't supported now!" << endl;
    }

    auto end = std::chrono::system_clock::now();
    auto endt_time = std::chrono::system_clock::to_time_t(end);

    // termination output
    if (is_normal)
        cout << "Normal termination of Jasp at "
             << ctime(&endt_time) << endl;
    else
        cout << "Error termination of Jasp at "
             << ctime(&endt_time) << endl;

    std::chrono::duration<double> elapsed_seconds = end-start;
    cout << "elapsed time: " << elapsed_seconds.count() << "s" << endl;


    return 0;
}
