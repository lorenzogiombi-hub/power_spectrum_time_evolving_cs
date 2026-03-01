// declare all functions

#ifndef FUNCTIONS_PGW_H
#define FUNCTIONS_PGW_H

// Declare the functions


std::vector<double> power_spectrum_mu_4D_interpolation(double kp, double eta_end, double omega, const std::vector<double>& k, const std::vector<double>& mu_vals, const std::vector<double>& x_vals, const std::vector<double>& tau1_vals, const std::vector<double>& tau2_vals, const std::vector<double>& x_data, const std::vector<double>& Pv_data);

#endif