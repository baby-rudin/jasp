#ifndef ENV_H
#define ENV_H

#include "type.hpp"
#include <string>
#include <vector>

extern std::string  CALC_METHOD;    // calculate method, default: RHF
extern std::string  BASIS_SETS;     // basis sets used, default: GEN

extern std::string  UNIT_LENGTH;    // unit of length, default: ANGSTROM
extern REAL         LEN_COEFF;      // change with LEN_UNIT
extern std::string  SCF_MIX;        // whether mix density last cycle, default: AUTO. MIX、NOMIX
extern REAL         SCF_E_EPS;      // delta E, default 1e-6
extern REAL         SCF_P_EPS;      // (sum( (P^{i}-P^{i-1})^2 ) / K^2 )^0.5, default: 1e-4
extern INTG         SCF_MAX_ROUND;  // max round of scf, default: 300

extern std::string  PRINT_LEVEL;    // control print what, default: ALL

// get command line  and  set environmental variables
bool  set_env(std::string line);

// output normal control information
void  print_ctrl_info();

#endif // ENV_H
