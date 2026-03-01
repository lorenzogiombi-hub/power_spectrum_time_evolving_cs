// declare all functions

#ifndef BESSEL_FUNCS_H
#define BESSEL_FUNCS_H

// Declare the functions
std::vector<double> linspace(double start, double end, int num);
std::vector<double> logspace(double start_exp, double end_exp, int num);

// double hypergeometric_0F1(double a, double z, int terms = 50);
double spherical_bessel_j(double nu, double x);
double spherical_bessel_y(double nu, double x);
double spherical_bessel_j_approx(double nu, double x);
double spherical_bessel_j_adaptive(double nu, double x);



#endif