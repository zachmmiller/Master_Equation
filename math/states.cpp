#include "states.h"

#include <math.h>

#include <algorithm>
#include <cstring>
#include <format>
#include <iostream>
#include <string>
#include <fstream>

#include "formulas.h"

Vibrational_Modes::Vibrational_Modes(int* vib_modes, int* vib_degen, double* ir_intens, double temperature, int N)
    : C(new int[N]), D(new int[N]), I(new double[N]), A(new double[N]), B(new double[N]), P(new double[N]), N(N) {
    memcpy((void*)this->C, (const void*)vib_modes, sizeof(int) * N);
    memcpy((void*)this->D, (const void*)vib_degen, sizeof(int) * N);
    memcpy((void*)this->I, (const void*)ir_intens, sizeof(double) * N);
    form::Compute_A_B_Coeffs(vib_modes, ir_intens, this->A, this->B, N, 1.0);
    form::Compute_Planck_Coeffs(vib_modes, this->P, N, temperature);
}

Vibrational_Modes::~Vibrational_Modes() {
    delete[] C;
    delete[] D;
    delete[] I;
    delete[] A;
    delete[] B;
    delete[] P;
}

void Vibrational_Modes::print() {
    std::cout << "Vibrational modes:\nIndex | mode (cm^-1) | degeneracy | IR intensity (km/mol) | A (s^-1) | B (cm^3 * J^-1 * s^-2) | B * P (s^-1)" << std::endl;
    for (int i = 0; i < N; i++) {
        std::cout << std::format("{:5d} | {:12d} | {:10d} | {:21.4f} | {:8.4f} | {:22.4e} | {:12.4f}", i, C[i], D[i], I[i], A[i], B[i], B[i] * P[i]) << std::endl;
    }
    std::cout << std::endl;
}

Occupation::Occupation(Vibrational_Modes* modes, long int energy) {
    N = modes->N;
    bounds = new int[N];
    int size = 0;
    for (long int i = 0; i < N; i++) {
        bounds[i] = (energy / modes->C[i]) + 1;
        size += bounds[i];
    }
    O = new double[size];
    std::fill(O, O + size, 0);
    this->energy = energy;
}

Occupation::~Occupation() {
    delete[] O;
    delete[] bounds;
}

int Occupation::offset(int index) {
    int out = 0;
    for (int i = 0; i < index; i++) {
        out += bounds[i];
    }
    return out;
}

double* Occupation::operator[](int index) {
    int offset = this->offset(index);
    return O + offset;
}

void Occupation::print() {
    double* mode;
    std::cout << "Density of states: " << density_of_states << std::endl;
    for (int i = 0; i < N; i++) {
        mode = (*this)[i];
        std::cout << "Mode " << i << " occupation:" << std::endl;
        for (int j = 0; j < bounds[i]; j++) {
            std::cout << "Level: " << j << ", occupation: " << mode[j] << ", probability: " << mode[j]/density_of_states << std::endl;
        }
    }
}

void Occupation::print(Vibrational_Modes* modes) {
    double* mode;
    std::cout << std::format("{:<20} {:15d}\n{:<20} {:15.8f}\n\n", "Energy (cm^-1):", energy, "Density of states:", density_of_states);
    std::string out_str;
    out_str += std::format("{:>20} | {:>20} | {:>20} | {:>20} | {:>20}\n", "Mode (cm^-1)", "Degeneracy", "Level", "Occupation", "Probability");
    for (int i = 0; i < N; i++) {
        mode = (*this)[i];
        for (int j = 0; j < bounds[i]; j++) {
            out_str += std::format("{:20d} | {:20d} | {:20d} | {:20.8f} | {:20.8f}\n", modes->C[i], modes->D[i], j, mode[j], mode[j]/density_of_states);
        }
        out_str += "\n";
        std::cout << out_str;
        out_str.clear();
    }
}

void Occupation::save(fs::path path, Vibrational_Modes* modes){
    double* mode;
    std::ofstream out_file;
    out_file.open(path);
    out_file << std::format("{:<20} {:15d}\n{:<20} {:15.8e}\n\n", "Energy (cm^-1):", energy, "Density of states:", density_of_states);
    std::string out_str;
    out_str += std::format("{:>20}, {:>20}, {:>20}, {:>20}, {:>20}\n", "Mode (cm^-1)", "Degeneracy", "Level", "Occupation", "Probability");
    for (int i = 0; i < N; i++) {
        mode = (*this)[i];
        for (int j = 0; j < bounds[i]; j++) {
            out_str += std::format("{:20d}, {:20d}, {:20d}, {:20.8e}, {:20.8e}\n", modes->C[i], modes->D[i], j, mode[j], mode[j]/density_of_states);
        }
        out_str += "\n";
        out_file << out_str;
        out_str.clear();
    }
}

void Occupation_Density(Vibrational_Modes* modes, Occupation* occupation, long int energy, double min_occupation) {
    energy++;
    double* density = new double[energy];
    std::fill(density, density + energy, 0);
    Beyer_Swinehart(modes->C, modes->D, modes->N, density, energy, -1);
    occupation->density_of_states = density[energy - 1];
    double* levels;
    double sum;

    for (int i = 0; i < modes->N; i++) {
        levels = (*occupation)[i];
        std::fill(density, density + energy, 0);
        Beyer_Swinehart(modes->C, modes->D, modes->N, density, energy, modes->C[i]);
        levels[1] = occupation->density_of_states - density[energy - 1];
        sum = 0;
        for (int j = 2; j < occupation->bounds[i]; j++) {
            std::fill(density, density + energy, 0);
            Beyer_Swinehart(modes->C, modes->D, modes->N, density, energy, modes->C[i] * j);
            levels[j] = occupation->density_of_states - density[energy - 1];
            levels[j - 1] -= levels[j];
            sum += levels[j - 1];
        }
        sum += levels[occupation->bounds[i] - 1];
        levels[0] = occupation->density_of_states - sum;
    }

    delete[] density;
}

void Beyer_Swinehart(int* C, int* D, int N, double* density, long int energy, long int eliminate_energy) {
    for (int i = 0; i < N; i++) {
        if (C[i] >= energy) {
            N = i;
            break;
        }
    }

    density[0] = 1;  // Because the ground state is populated

    if (eliminate_energy >= 0 && eliminate_energy < energy) {
        density[eliminate_energy]--;
    }

    int J, i, d, m;
    for (i = 0; i < N; i++) {
        J = C[i];
        for (d = 0; d < D[i]; d++) {
            density[J]++;
            for (m = J + 1; m < energy; m++) {
                density[m] = density[m] + density[m - J];
            }
        }
    }
}

double Density_of_States(int* C, int* D, int N, long int energy) {
    energy++;
    double* density = new double[energy];
    std::fill(density, density + energy, 0);
    Beyer_Swinehart(C, D, N, density, energy, -1);
    double density_of_states = density[energy - 1];
    delete[] density;
    return density_of_states;
}

double Sum_of_States(int* C, int* D, int N, long int energy) {
    energy++;
    double* density = new double[energy];
    std::fill(density, density + energy, 0);
    Beyer_Swinehart(C, D, N, density, energy, -1);
    double sum_of_states = 0;
    for (long int i = 0; i < energy; i++) {
        sum_of_states += density[i];
    }
    delete[] density;
    return sum_of_states;
}

void Boltzmann(std::complex<double>* vec, long int N_vec, int e_min, int e_step, int* v, int* D, int N_v, double T) {
    double* densities = new double[N_vec];
    double density_sum = 0;
    for (long int i = 0; i < N_vec; i++) {
        densities[i] = Density_of_States(v, D, N_v, e_step * i + e_min);
        density_sum += densities[i];
    }

    double* P = new double[N_vec];
    double P_sum = 0;
    for (long int i = 0; i < N_vec; i++) {
        P[i] = densities[i] * exp(-1 * (e_step * i + e_min) / (0.695034800 * T)) / density_sum;
        P_sum += P[i];
    }
    for (long int i = 0; i < N_vec; i++) {
        double& value_real = reinterpret_cast<double(&)[2]>(vec[i])[0];
        value_real = P[i] / P_sum;
    }

    delete[] densities;
    delete[] P;
}

/*
combo::combo(long int* data, long int N) {
    this->data = new long int[N];
    this->N = N;
    memcpy((void*)this->data, (const void*)data, sizeof(long int) * N);
}

combo::combo(combo& other) {
    this->data = new long int[other.N];
    this->N = other.N;
    memcpy((void*)this->data, (const void*)other.data, sizeof(long int) * other.N);
}

combo::combo(combo&& other) {
    this->data = other.data;
    this->N = other.N;
    other.data = nullptr;
    other.N = 0;
}

combo::~combo() { delete[] data; }

int combos::occurance(long int vibration, long int level) {
    int o = 0;
    for (combo& c : *this) {
        if (c.data[vibration] == level) {
            o++;
        }
    }
    return o;
}

long int combos::density(void) { return (long int)size(); }

Vibrations::Vibrations(long int* vibrations, long int* degeneracy, long int n_vibrations) : m_vibrations(vibrations), m_degeneracy(degeneracy), m_N_vibrations(n_vibrations) {
    m_degenerate_modes_size = 0;
    for (long int i = 0; i < n_vibrations; i++) {
        m_degenerate_modes_size += degeneracy[i];
    }
    m_degenerate_modes = new long int[m_degenerate_modes_size];

    long int offset = 0;
    for (long int i = 0; i < n_vibrations; i++) {
        for (long int j = 0; j < degeneracy[i]; j++) {
            m_degenerate_modes[offset] = vibrations[i];
            offset++;
        }
    }
}

Vibrations::~Vibrations() { delete[] m_degenerate_modes; }

long int Vibrations::offset(long int index) {
    long int out = 0;
    for (long int i = 0; i < index; i++) {
        out += m_degeneracy[i];
    }
    return out;
}
long int* Vibrations::get_degenerate_mode(long int index) {
    long int offset = this->offset(index);
    return m_degenerate_modes + offset;
}

void Vibrations::print() {
    std::cout << "Vibrations:";
    long int* mode;
    for (long int v = 0; v < m_N_vibrations; v++) {
        std::cout << std::endl;
        mode = get_degenerate_mode(v);
        for (long int s = 0; s < m_degeneracy[v]; s++) {
            std::cout << mode[s] << std::endl;
        }
    }
}

Occupations::Occupations(Vibrations* vibrations, long int energy) : m_vibrations(vibrations), m_energy(energy) {
    m_Ns = new long int[vibrations->m_N_vibrations];
    m_size = 0;
    for (long int i = 0; i < vibrations->m_N_vibrations; i++) {
        m_Ns[i] = (energy / vibrations->m_vibrations[i]) + 1;
        m_size += m_Ns[i];
    }
    m_occupations = new long int[m_size];
}
Occupations::~Occupations() {
    delete[] m_Ns;
    delete[] m_occupations;
}

long int Occupations::offset(long int index) {
    long int out = 0;
    for (long int i = 0; i < index; i++) {
        out += m_Ns[i];
    }
    return out;
}
long int* Occupations::operator[](long int index) {
    long int offset = this->offset(index);
    return m_occupations + offset;
}

long int Occupations::density(long int index) {
    long int* mode = (*this)[index];
    long int sum = 0;
    for (long int i = 0; i < m_Ns[index]; i++) {
        sum += mode[i];
    }
    return sum;
}

void Occupations::print() {
    long int* mode;
    long int density;
    for (long int v = 0; v < m_vibrations->m_N_vibrations; v++) {
        mode = (*this)[v];
        density = this->density(v);
        std::cout << "Vibrational energy: " << m_vibrations->m_vibrations[v] << ", degeneracy: " << m_vibrations->m_degeneracy[v] << std::endl;
        std::cout << "Density of states: " << density << ", scaled for degeneracy: " << density / m_vibrations->m_degeneracy[v] << std::endl;
        for (long int s = 0; s < m_Ns[v]; s++) {
            std::cout << "level: " << s << ", occupation: " << mode[s] << std::endl;
        }
    }
}



void Populate_Occupations(Vibrations* V, Occupations* O) {
    long int* S = new long int[V->m_N_vibrations];
    long int* mult = new long int[V->m_N_vibrations];
    long int M = O->m_Ns[0];
    long int K = V->m_N_vibrations;
    long int* C = V->m_vibrations;
    long int E = O->m_energy;
    long int* D = V->m_degeneracy;

    std::fill(S, S + K, 0);
    long int value;
    long int current_loop = K - 1;
    while (S[0] < M) {
        value = 0;
        for (long int i = 0; i < K; i++) {
            value += C[i] * S[i];
        }

        if (value < E) {
            current_loop = K - 1;
            S[current_loop]++;
            value += C[current_loop];
        }

        if (value == E) {
            long int total = 1;
            for (long int i = 0; i < K; i++) {
                if (D[i] > 1) {
                    mult[i] = Density_of_States(V->get_degenerate_mode(i), D[i], S[i] * C[i]);
                    total *= mult[i];
                } else {
                    mult[i] = 1;
                }
            }
            for (long int i = 0; i < K; i++) {
                if (D[i] > 1) {
                    Populate_Degenerate_Occupations(V->get_degenerate_mode(i), D[i], (*O)[i], S[i] * C[i], total / mult[i]);
                }
            }

            for (long int i = 0; i < K; i++) {
                if (D[i] == 1) {
                    (*O)[i][S[i]] += total;
                }
            }
        }

        if (value >= E) {
            if (current_loop > 0) {
                current_loop--;
            }
            for (long int i = current_loop + 1; i < K; i++) {
                S[i] = 0;
            }
            S[current_loop]++;
        }
    }
    delete[] S;
    delete[] mult;
}

void Populate_Degenerate_Occupations(long int* degenerate_vibrations, long int degeneracy, long int* occupations, long int E, long int scalar) {
    long int* S = new long int[degeneracy];
    long int M = (E / degenerate_vibrations[0]) + 1;

    std::fill(S, S + degeneracy, 0);
    long int value;
    long int current_loop = degeneracy - 1;
    while (S[0] < M) {
        value = 0;
        for (long int i = 0; i < degeneracy; i++) {
            value += degenerate_vibrations[i] * S[i];
        }

        if (value < E) {
            current_loop = degeneracy - 1;
            S[current_loop]++;
            value += degenerate_vibrations[current_loop];
        }

        if (value == E) {
            for (long int i = 0; i < degeneracy; i++) {
                occupations[S[i]] += scalar;
            }
        }

        if (value >= E) {
            if (current_loop > 0) {
                current_loop--;
            }
            for (long int i = current_loop + 1; i < degeneracy; i++) {
                S[i] = 0;
            }
            S[current_loop]++;
        }
    }
    delete[] S;
}

void Beyer_Swinehart(long int* C, long int K, long int*& P, long int N) {
    for (long int i = 0; i < K; i++) {
        if (C[i] >= N) {
            K = i;
            break;
        }
    }
    P = new long int[N];
    std::fill(P, P + N, 0);
    // P[0] = 1;
    long int J, i, m;
    for (i = 0; i < K; i++) {
        J = C[i];
        P[J] = P[J] + 1;
        for (m = J + 1; m < N; m++) {
            P[m] = P[m] + P[m - J];
        }

}
}

void Beyer_Swinehart_D(long int C, long int D, long int*& P, long int E) {
    long int N = (E / C) + 1;

    P = new long int[N];
    std::fill(P, P + N, 0);
    P[0] = 1;
    long int J, i, m;
    for (i = 0; i < D; i++) {
        J = 1;  // C[i] - reindex the array
        P[J] = P[J] + 1;
        for (m = J + 1; m < N; m++) {
            P[m] = P[m] + P[m - J];
        }
    }
}

void Beyer_Swinehart_P(long int* C, long int K, long int*& P, long int N, long int vibration, long int state) {
    for (long int i = 0; i < K; i++) {
        if (C[i] >= N) {
            K = i;
            break;
        }
    }
    // std::cout << "vibration: " << vibration << ", state_energy: " << state_energy << std::endl;

    P = new long int[N];
    std::fill(P, P + N, 0);
    // if (state == 1) {
    P[C[vibration] * state]--;
    //}
    // P[0] = 1;
    long int J, i, m;
    for (i = 0; i < K; i++) {
        J = C[i];
        P[J]++;
        for (m = J + 1; m < N; m++) {
            P[m] = P[m] + P[m - J];
        }
    }
}

void Beyer_Swinehart_degenerate(long int* C, long int* D, long int K, long int*& P, long int N, long int elim_N) {
    for (long int i = 0; i < K; i++) {
        if (C[i] >= N) {
            K = i;
            break;
        }
    }

    // std::cout << "vibration: " << vibration << ", state_energy: " << state_energy << std::endl;

    P = new long int[N];
    std::fill(P, P + N, 0);
    // if (state == 1) {
    if (elim_N >= 0 && elim_N < N) {
        P[elim_N]--;
    }
    //}
    // P[0] = 1;
    long int J, i, m, d;
    for (i = 0; i < K; i++) {
        J = C[i];
        for (d = 0; d < D[i]; d++) {
            P[J]++;
            for (m = J + 1; m < N; m++) {
                P[m] = P[m] + P[m - J];
            }
        }
    }
}

long int Sum_of_States(long int* C, long int K, long int E) {
    E++;
    long int* P;
    Beyer_Swinehart(C, K, P, E);
    long int sum = 0;
    for (int i = 0; i < E; i++) {
        sum += P[i];
    }
    delete[] P;
    return sum;
}

long int Density_of_States(long int* C, long int K, long int E) {
    E++;
    long int* P;
    Beyer_Swinehart(C, K, P, E);
    long int density = P[E - 1];
    delete[] P;
    return density;
}

long int Density_of_States_D(long int C, long int D, long int E) {
    long int* P;
    Beyer_Swinehart_D(C, D, P, E);
    long int density = P[(E / C)];
    delete[] P;
    return density;
}

long int Density_of_States_P(long int* C, long int K, long int E, long int vibration, long int state) {
    E++;
    long int* P;
    Beyer_Swinehart_P(C, K, P, E, vibration, state);
    long int density = P[E - 1];
    delete[] P;
    return density;
}

long int Density_of_States_degenerate(long int* C, long int* D, long int K, long int E, long int vibration, long int state) {
    E++;
    long int* P;
    Beyer_Swinehart_degenerate(C, D, K, P, E, -1);
    long int bucket_width = C[0];
    long int N_buckets = E / bucket_width + 1;
    long int density = P[N_buckets - 1];
    delete[] P;
    return density;
}

void combinations(combos& out, long int dE, long int* C, long int K) {
    long int density = Density_of_States(C, K, dE);
    out.reserve(density);

    long int* S = new long int[K];
    long int M = (dE / C[0]) + 1;

    std::fill(S, S + K, 0);
    long int V;
    long int current_loop = K - 1;
    long int combos_found = 0;

    while (S[0] < M) {
        V = 0;
        for (long int i = 0; i < K; i++) {
            V += C[i] * S[i];
        }

        if (V < dE) {
            current_loop = K - 1;
            S[current_loop]++;
            V += C[current_loop];
        }

        if (V == dE) {
            out.emplace_back(S, K);
            combos_found++;
        }

        if (V >= dE) {
            if (current_loop > 0) {
                current_loop--;
            }
            for (long int i = current_loop + 1; i < K; i++) {
                S[i] = 0;
            }
            S[current_loop]++;
        }
    }
    std::cout << "Number of combos found: " << combos_found << std::endl;
    delete[] S;
}

void combinations_d(combos& out, long int dE, long int* C, long int* D, long int K) {
    long int density = Density_of_States(C, K, dE);
    out.reserve(density);

    long int* S = new long int[K];
    long int M = (dE / C[0]) + 1;

    std::fill(S, S + K, 0);
    long int V;
    long int current_loop = K - 1;
    long int combos_found = 0;

    std::vector<long int*> degenerate_modes;
    for (long int i = 0; i < K; i++) {
        degenerate_modes.push_back(new long int[D[i]]);
        for (long int j = 0; j < D[i]; j++) {
            degenerate_modes[i][j] = C[i];
        }
    }

    while (S[0] < M) {
        V = 0;
        for (long int i = 0; i < K; i++) {
            V += C[i] * S[i];
        }

        if (V < dE) {
            current_loop = K - 1;
            S[current_loop]++;
            V += C[current_loop];
        }

        if (V == dE) {
            long int count = 1;
            for (long int i = 0; i < K; i++) {
                if (D[i] > 1) {
                    count *= Density_of_States(degenerate_modes[i], D[i], S[i] * C[i]);
                }
            }

            combos_found += count;
            for (long int i = 0; i < count; i++) {
                out.emplace_back(S, K);
                // combos_found++;
            }
        }

        if (V >= dE) {
            if (current_loop > 0) {
                current_loop--;
            }
            for (long int i = current_loop + 1; i < K; i++) {
                S[i] = 0;
            }
            S[current_loop]++;
        }
    }
    std::cout << "Number of combos found: " << combos_found << std::endl;
    delete[] S;
    for (long int* d : degenerate_modes) {
        delete[] d;
    }
}

void combos_new(long int* C, long int* D, long int K, long int E, Occupations& O) {
    long int* S = new long int[K];
    long int M = (E / C[0]) + 1;

    std::fill(S, S + K, 0);
    long int V;
    long int current_loop = K - 1;
    while (S[0] < M) {
        V = 0;
        for (long int i = 0; i < K; i++) {
            V += C[i] * S[i];
        }

        if (V < E) {
            current_loop = K - 1;
            S[current_loop]++;
            V += C[current_loop];
        }

        if (V == E) {
            // O.add(C, S, D, K);
        }

        if (V >= E) {
            if (current_loop > 0) {
                current_loop--;
            }
            for (long int i = current_loop + 1; i < K; i++) {
                S[i] = 0;
            }
            S[current_loop]++;
        }
    }
    delete[] S;
}

long int degenerate_combos(long* C, long int K, long E, long* out) {
    long int M = (E / C[0]) + 1;
    long int* S = new long int[K];
    std::fill(S, S + K, 0);
    long int V;
    long int current_loop = K - 1;

    long int permutations = 0;
    while (S[0] < M) {
        V = 0;
        for (long int i = 0; i < K; i++) {
            V += C[i] * S[i];
        }

        if (V < E) {
            current_loop = K - 1;
            S[current_loop]++;
            V += C[current_loop];
        }

        if (V == E) {
            permutations++;
            for (long int i = 0; i < K; i++) {
                out[S[i]]++;
            }
        }

        if (V >= E) {
            if (current_loop > 0) {
                current_loop--;
            }
            for (long int i = current_loop + 1; i < K; i++) {
                S[i] = 0;
            }
            S[current_loop]++;
        }
    }
    delete[] S;
    return permutations;
}

void Occupation(long int* C, double* O, long int K, long int E) {
    double sum = 0;
    int max[] = {11, 5, 2, 1};
    for (int c = 0; c < K; c++) {
        for (int i = 1; i <= max[c]; i++) {
            sum += exp(-1 * C[c] * i / E);
        }
    }
    std::cout << "Sum: " << sum << std::endl;
    for (int c = 0; c < K; c++) {
        std::cout << "Vibration: " << C[c] << std::endl;
        for (int i = 1; i <= max[c]; i++) {
            std::cout << "Occupation of i=" << i << ": " << exp(-1 * C[c] * i / E) << std::endl;
        }
    }
}

void Populate(double* M, long int* E, long int N, long int* C, double* U, long int K, double T) {
    // M : The matrix
    // E : Energy bins
    // N : M and E dimension
    // C : vibrational frequencies (in wavenumber)
    // U : transition dipoles corresponding to each vibration
    // T : Temperature (for Planck distribution)
    // The density of states for all energy bins

    double* O = new double[K];
    for (long int i = 0; i < N; i++) {
        // i is a column in the matrix.
        long int E_column = E[i];
        for (long int j = 0; j < N; j++) {
        }
    }
}

*/
