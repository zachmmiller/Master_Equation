#include "master_equation.h"

#include "formulas.h"
#include "states.h"

void J_Column_Blackbody(Vibrational_Modes* modes, int e_min, int e_step, double* matrix, long int col, long int N) {
    // long int energy_diff;
    double* element;
    double probability;

    Occupation occupations(modes, col * e_step + e_min);
    Occupation_Density(modes, &occupations, col * e_step + e_min);

    for (long int row = 0; row < N; row++) {
        element = matrix + col * N + row;
        double energy_diff = abs(col - row) * e_step;
        if (row > col) {
            // Under this condition, we are going up in energy
            for (int i = 0; i < modes->N; i++) {
                if (energy_diff - (double)e_step / 2 <= (double)modes->C[i] && (double)modes->C[i] < energy_diff + (double)e_step / 2) {
                    if (occupations.density_of_states == 0) {
                        // This means that everything is in the ground state.
                        *element += modes->P[i] * modes->B[i] * modes->D[i];
                    } else {
                        for (int j = 0; j < occupations.bounds[i]; j++) {
                            probability = occupations[i][j] / occupations.density_of_states;
                            *element += (modes->P[i] * modes->B[i]) * probability * (j + 1) * modes->D[i];
                        }
                    }
                }
            }
        }
        if (col > row) {
            // Under this condition, we are going down in energy
            for (int i = 0; i < modes->N; i++) {
                if (energy_diff - (double)e_step / 2 <= (double)modes->C[i] && (double)modes->C[i] < energy_diff + (double)e_step / 2) {
                    if (occupations.density_of_states == 0) {
                        // This means that everything is in the ground state.
                        // *element += modes->A[i] + modes->P[i] * modes->B[i] * modes->D[i];
                    } else {
                        for (int j = 0; j < occupations.bounds[i]; j++) {
                            probability = occupations[i][j] / occupations.density_of_states;
                            *element += (modes->A[i] + modes->P[i] * modes->B[i]) * probability * j * modes->D[i];
                            //*element += modes->A[i] * probability * j;
                        }
                    }
                }
            }
        }
    }

    element = matrix + col * N + col;

    for (long int row = 0; row < N; row++) {
        if (row == col) {
            continue;
        }
        *element -= matrix[col * N + row];
    }

    // RRKM rate constant could be computed here, but breaking them out into a separate function
}

void J_Column_RRKM(int* v_act, int* D_act, int N_act, int* v_trans, int* D_trans, int N_trans, int activation_energy, int e_min, int e_step, double* matrix, long int col, long int N) {
    int energy = e_step * col + e_min;
    if (energy > activation_energy) {
        double* element = matrix + N * col + col;
        double density_of_states_act = Density_of_States(v_act, D_act, N_act, (long int)energy);
        double sum_of_states_trans = Sum_of_States(v_trans, D_trans, N_trans, (long int)(energy - activation_energy));

        double k_RRKM = sum_of_states_trans / (density_of_states_act * form::h * (1.0 / (form::h * form::c * 100)));
        *element -= k_RRKM;
    }
}
