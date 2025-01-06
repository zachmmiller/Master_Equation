#ifndef FORMULAS_H
#define FORMULAS_H

#include <math.h>

#include <iostream>

namespace form {

const double h = 6.62607 * pow(10, -34);
const double e0 = 8.85419 * pow(10, -12);
const double c = 2.99792 * pow(10, 8);
const double kB = 1.38065 * pow(10, -23);
const double NA = 6.02214 * pow(10, 23);

void print_constants();

double Einstein_A(double v, double I);
double Einstein_B(double v, double I);
void Compute_A_B_Coeffs(int* v, double* I, double* A, double* B, int N, double s);

double Planck_Distribution(double T, double v);
void Compute_Planck_Coeffs(int* v, double* P, int N, double T);

}  // namespace form

#endif
