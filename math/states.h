#ifndef STATES_H
#define STATES_H

#include <complex>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>

namespace fs = std::filesystem;

struct Vibrational_Modes {
    Vibrational_Modes(int* vib_modes, int* vib_degen, double* ir_intens, double temperature, int N);
    ~Vibrational_Modes();
    void print();
    int* C;
    int* D;
    double* I;  // IR Intensity
    double* A;  // Einstein A coefficients
    double* B;  // Einstein B coefficients
    double* P;  // Planck distribution
    int N;
};

struct Occupation {
    Occupation(Vibrational_Modes* modes, long int energy);
    ~Occupation();
    int offset(int index);
    double* operator[](int index);
    void print();
    void print(Vibrational_Modes* modes);
    void save(fs::path path, Vibrational_Modes* modes);
    double density_of_states;
    double* O;
    int* bounds;
    int N;
    long int energy;
};

double Planck_Distribution(double wavenumber, double temperature);

void Occupation_Density(Vibrational_Modes* modes, Occupation* occupation, long int energy, double min_occupation=0);
void Beyer_Swinehart(int* C, int* D, int N, double* density, long int energy, long int eliminate_energy);

double Density_of_States(int* C, int* D, int N, long int energy);
double Sum_of_States(int* C, int* D, int N, long int energy);

void Boltzmann(std::complex<double>* vec, long int N_vec, int e_min, int e_step, int* v, int* D, int N_v, double T);

/*
struct combo {
    combo(long int* data, long int N);
    combo(combo& other);
    combo(combo&& other);
    ~combo();
    long int* data;
    long int N;
};

class combos : public std::vector<combo> {
   public:
    int occurance(long int i, long int level);
    long int density(void);
};

struct Vibrations {
    Vibrations(long int* vibrations, long int* degeneracy, long int n_vibrations);
    ~Vibrations();
    long int* m_vibrations;
    long int* m_degeneracy;
    long int m_N_vibrations;
    long int* m_degenerate_modes;
    long int m_degenerate_modes_size;

    long int* get_degenerate_mode(long int index);
    long int offset(long int index);
    void print();
};

struct Occupations {
    Occupations(Vibrations* vibrations, long int energy);
    ~Occupations();

    long int* operator[](long int index);
    long int density(long int index);
    void print();

    long int offset(long int index);

    Vibrations* m_vibrations;
    long int m_energy;
    long int* m_Ns;
    long int m_size;
    long int* m_occupations;
};

void Beyer_Swinehart(long int* C, long int K, long int*& P, long int N);
void Beyer_Swinehart_D(long int C, long int D, long int*& P, long int E);
void Beyer_Swinehart_P(long int* C, long int K, long int*& P, long int N, long int vibration, long int state);
void Beyer_Swinehart_degenerate(long int* C, long int* D, long int K, long int*& P, long int N, long int elim_N);
long int Sum_of_States(long int* C, long int K, long int E);
long int Density_of_States(long int* C, long int K, long int E);
long int Density_of_States_D(long int C, long int D, long int E);
long int Density_of_States_P(long int* C, long int K, long int E, long int vibration, long int state);
long int Density_of_States_degenerate(long int* C, long int* D, long int K, long int E, long int vibration, long int state);

void Populate_Occupations(Vibrations* V, Occupations* O);
void Populate_Degenerate_Occupations(long int* degenerate_vibrations, long int degeneracy, long int* occupations, long int E, long int scalar);

void combinations(combos& out, long int dE, long int* C, long int K);
void combinations_d(combos& out, long int dE, long int* C, long int* D, long int K);

// void combos_new(long int* C, long int* D, long int K, long int E, Occupations& O);
long int degenerate_combos(long int* C, long int K, long int E, long int* out);

// void Occupation(long int* C, double* O, long int K, long int E);
void Populate(double* M, long int* E, long int N, long int* C, double* U, double T);
*/

#endif
