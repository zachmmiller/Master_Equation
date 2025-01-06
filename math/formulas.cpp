#include "formulas.h"

namespace form {

void print_constants() {
    std::cout << "Constants and units: " << std::endl;
    std::cout << "Planck constant (h): " << h << " (J * hz^-1)" << std::endl;
    std::cout << "Epsilon nought (e0): " << e0 << " (C^2 * kg^-1 * m^-3 * s^2)" << std::endl;
    std::cout << "Speed of light in vacuum (c): " << c << " (m * s^-1)" << std::endl;
    std::cout << "Boltzmann constant: " << kB << " (J * K^-1)" << std::endl;
}

double Einstein_A(double v, double I, double s) {
    // Input units:
    //    v: cm^-1
    //    I: KM / mol
    //    s: no units (this is a scaling factor for IR intensities)
    // Output units:
    //    s^-1

    v = v * 100;       // convert cm^-1 to m^-1
    I = I * 1000 * s;  // convert km/mol to m/mol and scale by s

    return 8 * M_PI * c * pow(v, 2) * I / NA;
}

double Einstein_B(double v, double A) {
    // Input units:
    //    v: cm^-1
    //    A: s^-1
    // Output units:
    //    cm^3 * s^-2 * J^-1

    double c_cm = c * 100;  // convert c from m/s to cm/s
    v = v * c_cm;           // convert v from cm^-1 to Hz

    // Only using cm here because that is what I have seen in the literature.
    // We end up with a rate that has units of s^-1 when we multiply B with the planck distribution,
    // so it really only matters that we are consistent here with Planck_Distribution below.

    return A * pow(c_cm, 3) / (8 * M_PI * h * pow(v, 3));
}

void Compute_A_B_Coeffs(int* v, double* I, double* A, double* B, int N, double s) {
    for (int i = 0; i < N; i++) {
        A[i] = Einstein_A(v[i], I[i], s);
        B[i] = Einstein_B((double)v[i], A[i]);
    }
}

double Planck_Distribution(double T, double v) {
    // Input units:
    //    T: k
    //    v: cm^-1
    // Output units:
    //    J * s * cm^-3

    double c_cm = c * 100;  // convert c from m/s to cm/s
    v = v * c_cm;           // convert v from cm^-1 to Hz

    return 8 * M_PI * h * pow(v, 3) / pow(c_cm, 3) * (1.0 / (exp(h * v / (kB * T)) - 1));
}

void Compute_Planck_Coeffs(int* v, double* P, int N, double T) {
    for (int i = 0; i < N; i++) {
        P[i] = Planck_Distribution(T, (double)v[i]);
    }
}

}  // namespace form
