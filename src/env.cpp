#include "env.h"
#include "stropt.h"
#include "constant.h"
#include <algorithm>
#include <iostream>

using namespace std;

string      CALC_METHOD    ("RHF");     // calculate method, default: RHF
string      BASIS_SETS     ("GEN");     // basis sets used, default: GEN
string      UNIT_LENGTH    ("ANG");     // unit of length , default: ANGSTROM
string      SCF_MIX        ("AUTO");    // whether mix density last cycle, default: auto
REAL        SCF_E_EPS      = 1e-8;      // delta E, default 1e-10
REAL        SCF_P_EPS      = 1e-6;      // (sum( (P^{i}-P^{i-1})^2 ) / K^2 )^0.5, default: 1e-8
INTG        SCF_MAX_ROUND  = 1000;      // max round of scf, default: 300

string      PRINT_LEVEL     ("NORM");   // control print what ALL, MORE, NORM, LESS

REAL        LEN_COEFF      = ANGSTROM_PER_BOHR;


bool    set_env(string line)
{
    bool    ret = true;

    replace( line.begin(), line.end(), '#', ' ');
    replace( line.begin(), line.end(), '/', ' ');
    clean_line(line);
    str_upper(line);

    auto commands = split_string(line, string(" "));

    CALC_METHOD = commands[0];
    BASIS_SETS  = commands[1];

    for (size_t i=2; i<commands.size(); i++) {
        if (commands[i].empty()) continue;

        auto cmd = split_string(commands[i], string("="));
        if (cmd[0] == "UNIT") {
            UNIT_LENGTH = cmd[1];
            if (UNIT_LENGTH == "BOHR")
                LEN_COEFF = REAL(1.0);
        }
        else if (cmd[0] == "SCF") {
            if (cmd.size() > 1)
                SCF_MIX = cmd[1];
        }
        else
        {
            ret = false;
        }
    }

    if (!ret) {
        cout << "Error control line!" << endl;
        exit(127);
    }

    return ret;
}


void  print_ctrl_info()
{

    string L1 =  string("   METHOD       = ") + CALC_METHOD;
    string L2 =  string("   BASIS        = ") + BASIS_SETS;
    string L3 =  string("   UNIT         = ") + UNIT_LENGTH ;

    cout << "      Route Information:" << endl;
    cout << " " + string(30, '-') << endl;
    cout << " " + L1 << endl;
    cout << " " + L2 << endl;
    cout << " " + L3 << endl;
    cout << " " + string(30, '-') << endl;
    cout << endl;
}

