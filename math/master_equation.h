#ifndef MASTER_EQUATION_H
#define MASTER_EQUATION_H

#include "states.h"

void J_Column_Blackbody(Vibrational_Modes* modes, int e_min, int e_step, double* matrix, long int col, long int N);

void J_Column_RRKM(int* v_act, int* D_act, int N_act, int* v_trans, int* D_trans, int N_trans, int activation_energy, int e_min, int e_step, double* matrix, long int col, long int N);

#endif
